#include <cstdint>
#include <iostream>
#include <cassert>
#include <cstring>
#include <sstream>
#include <stdexcept>
#include <cstdint>
#include <math.h>

#include <algorithm>
#include <map>
#include <vector>

#include "kff_io.hpp"

using namespace std;


// Utils

template<typename T>
void store_big_endian(uint8_t * buff, size_t size, const T& data) {
	uint64_t value = data;
	for (int b = size - 1; b >= 0; --b) {
		*buff++ = value >> (8 * b);
	}
}

template<typename T>
void load_big_endian(uint8_t * buff, size_t size, T& data) {
	data = 0;
	for (uint b=0 ; b < size; b++) {
		data <<= 8;
		data |= buff[b];
	}
}

uint64_t bytes_from_bit_array(uint64_t bits_per_elem, uint64_t nb_elem) {
	if (bits_per_elem == 0 or nb_elem == 0)
		return 0;
	else
		return ((bits_per_elem * nb_elem - 1) / 8) + 1;
}

static void rightshift8(uint8_t * bitarray, size_t length, size_t bitshift);


/* Copy nb nucleotides of src, starting at nucleotide src_pos, into dst at nucleotide dst_pos.
 * The nucleotide i of an array is stored in the 2 bits 7-2(i%4) and 6-2(i%4) of its byte i/4
 * (big endian order). The bits of the copied area must be 0 in dst (bitwise or).
 */
static void copy_nucleotides(uint8_t * dst, uint64_t dst_pos, const uint8_t * src, uint64_t src_pos, uint64_t nb) {
	if (nb == 0)
		return;
	src += src_pos / 4;
	dst += dst_pos / 4;
	const unsigned src_bit = 2 * (src_pos % 4);
	const unsigned dst_bit = 2 * (dst_pos % 4);
	const uint64_t bits = 2 * nb;
	const uint64_t src_bytes = (src_bit + bits + 7) / 8;
	const uint64_t dst_bytes = (dst_bit + bits + 7) / 8;
	// Mask of the bits of the last destination byte that are part of the copy
	const unsigned last_bits = (dst_bit + bits) - 8 * (dst_bytes - 1);
	const uint8_t last_mask = (uint8_t)(0xFF << (8 - last_bits));

	if (src_bit == dst_bit) {
		// Same alignment: whole bytes are copied
		if (dst_bytes == 1) {
			dst[0] |= src[0] & (0xFF >> dst_bit) & last_mask;
			return;
		}
		dst[0] |= src[0] & (0xFF >> dst_bit);
		memcpy(dst + 1, src + 1, dst_bytes - 2);
		dst[dst_bytes - 1] |= src[dst_bytes - 1] & last_mask;
		return;
	}

	// 8 bits of the source from the bit p (only the bytes of the source are read)
	auto get8 = [src, src_bytes](uint64_t p) -> uint8_t {
		const uint64_t byte = p / 8;
		const unsigned offset = p % 8;
		unsigned word = (unsigned)src[byte] << 8;
		if (offset != 0 and byte + 1 < src_bytes)
			word |= src[byte + 1];
		return (uint8_t)(word >> (8 - offset));
	};

	// First destination byte: bits dst_bit..7 receive the source bits from src_bit
	uint8_t first = dst_bit > src_bit ? (uint8_t)(get8(0) >> (dst_bit - src_bit)) : get8(src_bit - dst_bit);
	uint8_t first_mask = 0xFF >> dst_bit;
	if (dst_bytes == 1)
		first_mask &= last_mask;
	dst[0] |= first & first_mask;
	if (dst_bytes == 1)
		return;
	// Following bytes: the destination bit 8j comes from the source bit src_bit + 8j - dst_bit
	for (uint64_t j = 1; j + 1 < dst_bytes; j++)
		dst[j] = get8(src_bit + 8 * j - dst_bit);
	dst[dst_bytes - 1] |= get8(src_bit + 8 * (dst_bytes - 1) - dst_bit) & last_mask;
}


// ----- Open / Close functions -----

Kff_file::Kff_file(const string filename, const string mode) {
	// Variable init
	this->filename = filename;
	
	this->is_writer = false;
	this->is_reader = false;
	
	this->writing_started = false;
	this->next_free = 0;
	this->buffer_size = 1 << 10; // 1 KB
	// this->buffer_size = 1 << 4;
	this->max_buffer_size = 1 << 20; // 1 MB
	// this->max_buffer_size = 1 << 6;
	this->file_buffer = new uint8_t[this->buffer_size];
	this->file_size = 0;
	this->delete_on_destruction = false;
	this->footer = nullptr;

	try {
		this->open(mode);
	} catch (...) {
		// The destructor is not called when the constructor throws
		delete[] this->file_buffer;
		delete this->footer;
		for (Section_Index * si : this->index)
			delete si;
		throw;
	}
}

void Kff_file::open(string mode) {
	this->writing_started = false;
	this->current_position = 0;
	this->read_cache_start = 0;
	this->read_cache_size = 0;

	// Determine the mode and open the file
	if (mode[0] == 'w') {
		this->is_writer = true;
		this->file_size = 0;
		this->next_free = 0;
	} else if (mode[0] == 'r') {
		this->is_reader = true;
		// If no info on the file
		if (this->file_size == 0 and this->next_free == 0) {
			// Open the fp
			this->fs.open(this->filename, fstream::binary | fstream::in);
			if (this->fs.fail()) {
				std::string msg = "Cannot open file " + this->filename;
				throw std::runtime_error(msg);
			}
			// Compute the file length
			long position = this->fs.tellp();
			this->fs.seekg(0, this->fs.end);
			this->file_size = (long)(this->fs.tellp()) - position;
			// Go back to the beginning
			this->fs.seekg(0, this->fs.beg);
		}
	} else {
		cerr << "Unsupported mode " << mode << endl;
		exit(1);
	}

	this->tmp_closed = false;
	this->header_over = false;
	this->indexed = false;
	this->footer = nullptr;
	this->footer_discovery_ended = true;

	// Write the signature and the version at the beginning of the file
	if (this->is_writer) {
		uint8_t default_encoding = 0b00011110;
		// Signature
		uint8_t buff[] = {	'K', 'F', 'F',
							KFF_VERSION_MAJOR, KFF_VERSION_MINOR,
							default_encoding,
							0 /*uniqueness*/, 0 /*canonicity*/
						};

		this->write(buff, 8);

		this->indexed = true;
		this->end_position = 0;
	}
	// Read the header
	else if (this->is_reader) {
		// Header integrity marker
		uint8_t buff[4];
		this->read(buff, 3);
		if (buff[0] != 'K' or buff[1] != 'F' or buff[2] != 'F') {
			cerr << "Absent KFF signature at the beginning of the file." << endl;
			cerr << "Please check that the file is not corrupted" << endl;
			throw "Absent signature at the beginning";
		}

		// Version reading
		this->read(&this->major_version, 1);
		this->read(&this->minor_version, 1);
		if (KFF_VERSION_MAJOR < this->major_version or (KFF_VERSION_MAJOR == this->major_version and KFF_VERSION_MINOR < this->minor_version)) {
			cerr << "The software version " << (uint)KFF_VERSION_MAJOR << "." << (uint)KFF_VERSION_MINOR << " can't read files writen in version " << (uint)this->major_version << "." << (uint)this->minor_version << endl;
			throw "Unexpected version number";
		}
		// Encoding load
		this->read_encoding();
		// Read global flags
		uint8_t flag;
		this->read(&flag, 1);
		this->uniqueness = flag != 0;

		this->read(&flag, 1);
		this->canonicity = flag != 0;
		// Read metadata size
		this->read(buff, 4);
		load_big_endian(buff, 4, this->metadata_size);


		// Footer integrity marker
		unsigned long saved_position = this->tellp();
		this->jump_to(3, true);
		this->end_position = this->tellp();
		this->read(buff, 3);
		this->jump_to(saved_position);
		if (buff[0] != 'K' or buff[1] != 'F' or buff[2] != 'F') {
			cerr << "Absent KFF signature at the end of the file." << endl;
			cerr << "Please check that the file is not corrupted" << endl;
			throw "Absent signature at the end";
		}

		// Back to the start
		this->footer_discovery_ended = false;
		// Discover footer
		this->footer_discovery();
		this->index_discovery();
	}
}

void Kff_file::close(bool write_buffer) {
	if (this->is_writer) {
		// Write the index
		if (this->indexed)
			this->write_footer();
		// Write the signature
		char signature[] = {'K', 'F', 'F'};
		this->write((uint8_t *)signature, 3);
		
		// Write the end of the file
		if (write_buffer) {
			// The file was never opened
			if (not this->writing_started) {
				this->writing_started = true;
				this->fs.open(this->filename, fstream::binary | fstream::out);
			} else if (this->tmp_closed) {
				this->reopen();
			}
			// Write the buffer
			this->fs.write((char *)this->file_buffer, this->next_free);
			if (this->fs.fail()) {
				cerr << "Filesystem problem during buffer disk saving" << endl;
				exit(1);
			}
			this->file_size += this->next_free;
			this->next_free = 0;
		} else {
			this->delete_on_destruction = true;
		}

		// cout << "delete_on_destruction " << delete_on_destruction << endl;
		// cout << this->filename << endl;
		if (this->fs.is_open())
			this->fs.close();
	}
	else if (this->is_reader) {
		// Release the read cache (the object can stay alive a long time, e.g. in merge)
		std::vector<uint8_t>().swap(this->read_cache);
		this->read_cache_start = 0;
		this->read_cache_size = 0;
	}

	this->tmp_closed = false;
	this->is_writer = false;
	this->is_reader = false;
}


Kff_file::~Kff_file() {
	this->close();

	delete[] this->file_buffer;
	if (this->delete_on_destruction and this->file_size > 0) {
		remove(this->filename.c_str());
	}

	if (this->footer != nullptr)
		delete this->footer;

	for (Section_Index * si : this->index)
		delete si;
}


void Kff_file::set_indexation(bool indexed) {
	if (this->is_writer)
		this->indexed = indexed;
}


void Kff_file::register_position(char section_type) {
	if (this->is_writer and this->indexed) {
		this->section_positions[this->tellp()] = section_type;
	}
}


void Kff_file::complete_header() {
	if (this->header_over)
		return;

	// If the metadata has not been read, jump over
	if (this->is_reader) {
		this->jump(this->metadata_size);
	}

	// If metadata has not been write, write a 0 byte one.
	else if (this->is_writer) {
		this->write_metadata(0, nullptr);
	}

	this->header_over = true;
}


void Kff_file::footer_discovery() {
	long current_pos = this->tellp();

	// Look at the footer: the file ends with "footer_size\0" + 8 bytes + "KFF".
	// The footer is optional, small files cannot contain one.
	bool footer_present = false;
	if (this->file_size + this->next_free >= 23) {
		this->jump_to(23, true);
		char name[11];
		this->read((uint8_t *)name, 11);
		footer_present = memcmp(name, "footer_size", 11) == 0;
	}

	if (footer_present) {
		this->jump(1); // remove the '\0'

		uint64_t size = 0;
		uint8_t buff[8];
		this->read(buff, 8);
		load_big_endian(buff, 8, size);
		// Jump to value section start
		this->jump_to(size+3, true);

		this->footer = new Section_GV(this);
		this->footer->close();
	}

	this->footer_discovery_ended = true;
	this->jump_to(current_pos);
}


void Kff_file::index_discovery() {
	long current_pos = this->tellp();
	bool header_over = this->header_over;
	this->complete_header();

	// Search in footer
	if (this->footer != nullptr and this->footer->vars.find("first_index") != this->footer->vars.end()) {
		this->indexed = true;
		this->read_index((long)this->footer->vars["first_index"]);
	}

	// Search first section
	if (not this->indexed) {
		char type = (char)this->peek_byte();
		if (type == 'i') {
			this->indexed = true;
			this->read_index(this->tellp());
		}

	}

	this->header_over = header_over;
	this->index_discovery_ended = true;

	this->jump_to(current_pos);
}


void Kff_file::read_index(long position) {
	long init_pos = this->tellp();

	while (position != 0) {
		// Move to the beginning
		this->jump_to(position);
		// read the local index content
		Section_Index * si = new Section_Index(this);
		this->index.push_back(si);
		si->close();
		// Update index position to the next index section
		if (si->next_index == 0)
			position = 0;
		else {
			position = this->tellp() + si->next_index;
		}
	}

	this->jump_to(init_pos);
}


// Size of the read cache. Moderate because a reader is created for each input file of a merge.
static const unsigned long READ_CACHE_CAPACITY = 1 << 16; // 64 KB

void Kff_file::read(uint8_t * bytes, unsigned long size) {
	if (not this->is_reader) {
		cerr << "Cannot read a file in writing mode." << endl;
		exit(1);
	}

	// Fast path: all the bytes are in the read cache
	if (this->current_position >= this->read_cache_start and
		this->current_position + size <= this->read_cache_start + this->read_cache_size) {
		if (size > 0)
			memcpy(bytes, this->read_cache.data() + (this->current_position - this->read_cache_start), size);
		this->current_position += size;
		return;
	}

	this->read_uncached(bytes, size);
}

void Kff_file::read_uncached(uint8_t * bytes, unsigned long size) {
	// Part of the bytes stored in the file on disk
	while (size > 0 and this->current_position < this->file_size) {
		const unsigned long cache_end = this->read_cache_start + this->read_cache_size;
		if (this->current_position >= this->read_cache_start and this->current_position < cache_end) {
			// Copy what is available in the cache
			unsigned long nb = std::min(size, cache_end - this->current_position);
			memcpy(bytes, this->read_cache.data() + (this->current_position - this->read_cache_start), nb);
			bytes += nb;
			size -= nb;
			this->current_position += nb;
		} else if (size >= READ_CACHE_CAPACITY) {
			// Large reads bypass the cache
			unsigned long nb = std::min(size, this->file_size - this->current_position);
			this->disk_read(bytes, nb, this->current_position);
			bytes += nb;
			size -= nb;
			this->current_position += nb;
		} else {
			this->fill_read_cache(this->current_position);
		}
	}

	// Part of the bytes in the memory buffer (file not entirely written on disk)
	if (this->current_position >= this->file_size) {
		uint64_t buffer_position = this->current_position - this->file_size;
		if (buffer_position + size > this->next_free) {
			string error = string("Read out of the file, Byte ") + to_string(this->file_size + this->next_free);
			throw out_of_range(error);
		}

		if (size > 0)
			memcpy(bytes, this->file_buffer + buffer_position, size);
		this->current_position += size;
	}
}

void Kff_file::disk_read(uint8_t * bytes, unsigned long size, unsigned long position) {
	// File not opened
	if (not this->fs.is_open())
		this->fs.open(this->filename, fstream::binary | fstream::in);

	this->fs.seekg(position);
	this->fs.read((char *)bytes, size);
	if (this->fs.fail()) {
		cerr << "Impossible to read the file " << this->filename << " on disk." << endl;
		exit(1);
	}
}

void Kff_file::fill_read_cache(unsigned long position) {
	if (this->read_cache.size() < READ_CACHE_CAPACITY)
		this->read_cache.resize(READ_CACHE_CAPACITY);

	unsigned long nb = std::min(READ_CACHE_CAPACITY, this->file_size - position);
	// Invalidate the cache before reading (disk_read can stop the program)
	this->read_cache_size = 0;
	this->disk_read(this->read_cache.data(), nb, position);
	this->read_cache_start = position;
	this->read_cache_size = nb;
}

uint8_t Kff_file::peek_byte() {
	// Byte on disk
	if (this->current_position < this->file_size) {
		if (this->current_position < this->read_cache_start or
			this->current_position >= this->read_cache_start + this->read_cache_size)
			this->fill_read_cache(this->current_position);
		return this->read_cache[this->current_position - this->read_cache_start];
	}
	// Byte in the memory buffer
	if (this->current_position < this->file_size + this->next_free)
		return this->file_buffer[this->current_position - this->file_size];
	// End of the file
	return 0;
}

void Kff_file::write(const uint8_t * bytes, unsigned long size) {
	if (not this->is_writer) {
		if (this->is_reader)
			cerr << "Cannot write a file in reading mode." << endl;
		else
			cerr << "Cannot write a closed file" << endl;
		exit(1);
	}

	unsigned long buff_space = this->buffer_size - this->next_free;

	// Resize buffer
	while (buff_space < size and this->buffer_size < this->max_buffer_size) {
		// Enlarge the buffer
		this->buffer_size *= 2;
		uint8_t * next_buffer = new uint8_t[this->buffer_size];
		// Copy the previous values
		memcpy(next_buffer, this->file_buffer, this->next_free);
		buff_space = this->buffer_size - this->next_free;
		// Fill the empty part with 0s
		memset(next_buffer + this->next_free, 0, buff_space);
		// remove the previous space
		delete[] this->file_buffer;
		this->file_buffer = next_buffer;
	}

	// fill the buffer
	if (buff_space >= size) {
		if (size > 0)
			memcpy(this->file_buffer + this->next_free, bytes, size);
		this->next_free += size;
	}
	// Not enought space, write the file
	else {
		// Open the file if needed
		if (not this->writing_started) {
			this->fs.open(this->filename, fstream::binary | fstream::out);
			this->writing_started = true;
		} else if (this->tmp_closed) {
			this->reopen();
		}

		this->fs.write((char*)this->file_buffer, this->next_free);
		this->fs.write((char*)bytes, size);
		this->file_size += this->next_free + size;
		this->next_free = 0;

		if (this->fs.fail()) {
			cerr << "File system error while writing " << this->filename << endl;
			exit(1);
		}
	}

	this->current_position += size;
}

void Kff_file::write_at(const uint8_t * bytes, unsigned long size, unsigned long position) {
	if (not this->is_writer) {
		if (this->is_reader)
			cerr << "Cannot write a file in reading mode." << endl;
		else
			cerr << "Cannot write a closed file" << endl;
		exit(1);
	}

	if (position > this->file_size + this->next_free) {
		cerr << "Cannot write after the last byte of the file." << endl;
		exit(1);
	}

	// Write the file on disk
	if (position < this->file_size) {
		// Only in file
		if (position + size <= this->file_size) {
			if (this->tmp_closed) {
				this->reopen();
			}
			this->fs.seekp(position);
			this->fs.write((char*)bytes, size);
			if (this->fs.fail()) {
				cerr << "File system error while writing " << this->filename << " at position " << position << endl;
				exit(1);
			}
			this->fs.seekp(this->file_size);
		}
		// On both file and buffer
		else {
			unsigned long in_file_size = this->file_size - position;
			// Write the file part
			this->write_at(bytes, in_file_size, position);
			// Write the buffer part
			this->write_at(bytes + in_file_size, size - in_file_size, position + in_file_size);
		}
	}
	// Write the buffer in RAM
	else {
		unsigned long corrected_position = position - this->file_size;
		
		// Write in the current buffer space
		if (corrected_position + size <= this->next_free) {
			memcpy(this->file_buffer + corrected_position, bytes, size);
		}
		// Spillover the buffer
		else {
			this->next_free = corrected_position;
			this->write(bytes, size);
		}
	}
}

unsigned long Kff_file::tellp() {
	return this->current_position;
}


void Kff_file::jump(long size) {
	// cout << "Jump " << this->current_position << " " << size << " / " << this->file_size << " " << this->next_free << endl;
	this->jump_to(this->current_position + size);
}

void Kff_file::jump_to(unsigned long position, bool from_end) {
	if (this->file_size + this->next_free < position) {
		cerr << "Jump out of the file." << endl;
		exit(1);
	}

	// Determine absolute position
	if (from_end) {
		position = this->file_size + this->next_free - position;
	}

	// In reading mode, the read cache is used: the stream is moved only when the cache is refilled
	if (this->is_reader) {
		this->current_position = position;
		return;
	}

	// Jump into the written file
	if (position < this->file_size) {
		this->fs.seekp(position);
	}
	// Jump into the buffer
	else /*if (this->current_position < this->file_size)*/ {
		this->fs.seekg(0, this->fs.end);
	}
	this->current_position = position;
}


void Kff_file::tmp_close() {
	if (this->is_writer and this->fs.is_open()) {
		this->fs.close();
		this->fs.clear();
		this->tmp_closed = true;
	}
}


void Kff_file::reopen() {
	if (this->tmp_closed) {
		auto streammode = fstream::binary | fstream::out | fstream::in | fstream::ate;

		// Open the file
		this->fs.open(this->filename, streammode);
		this->tmp_closed = false;
	}
}


void Kff_file::write_footer() {
	Section_Index si(this);

	// Compute end position
	long position = si.beginning + 17 + 9 * this->section_positions.size();
	// Add the values
	for (auto & it : this->section_positions) {
		si.register_section(it.second, it.first - position);
	}

	si.close();

	// Write a value section to register everything
	Section_GV sgv(this);
	sgv.write_var("first_index", si.beginning);
	sgv.write_var("footer_size", 9 + 2 * (12 + 8));
	sgv.close();
}


// ----- Header functions -----

void Kff_file::write_encoding(uint8_t a, uint8_t c, uint8_t g, uint8_t t) {
	// Value masking
	a &= 0b11; c &= 0b11; g &= 0b11; t &= 0b11;

	// Verification of the differences (explicit check: must not depend on NDEBUG).
	if (a == c or a == g or a == t or c == g or c == t or g == t)
		throw std::invalid_argument("Wrong encoding. The 4 2-bits values must be different.");

	// set values
	this->encoding[0] = a;
	this->encoding[1] = c;
	this->encoding[2] = g;
	this->encoding[3] = t;

	// Write to file
	uint8_t code = (a << 6) | (c << 4) | (g << 2) | t;
	this->write_at(&code, 1, 5);
}

void Kff_file::set_uniqueness(bool uniqueness) {
	uint8_t bit_uniq = uniqueness ? 1 : 0;
	this->write_at(&bit_uniq, 1, 6);
}
void Kff_file::set_canonicity(bool canonicity) {
	uint8_t bit_canon = canonicity ? 1 : 0;
	this->write_at(&bit_canon, 1, 7);
}

void Kff_file::write_encoding(uint8_t * encoding) {
	this->write_encoding(encoding[0], encoding[1], encoding[2], encoding[3]);
}

void Kff_file::read_encoding() {
	uint8_t code, a, c, g, t;
	// Get code values
	this->read(&code, 1);

	// Split each nucleotide encoding
	this->encoding[0] = a = (code >> 6) & 0b11;
	this->encoding[1] = c = (code >> 4) & 0b11;
	this->encoding[2] = g = (code >> 2) & 0b11;
	this->encoding[3] = t = code & 0b11;

	// Verification of the encoding
	if (a == c or a == g or a == t or c == g or c == t or g == t) {
		throw "Wrong encoding. The 4 2-bits values must be different.";
	}
}

void Kff_file::write_metadata(uint32_t size, const uint8_t * data) {
	if (this->header_over) {
		cerr << "The metadata have to be written prior to other content." << endl;
		exit(1);
	}

	uint8_t buff[4];
	store_big_endian(buff, 4, size);
	this->write(buff, 4);
	this->write(data, size);

	this->header_over = true;
}


void Kff_file::read_metadata(uint8_t * data) {
	this->read(data, this->metadata_size);
	this->header_over = true;
}

bool Kff_file::jump_next_section() {
	if (not is_reader)
		return false;
	char section_type = read_section_type();
	if (this->current_position == this->file_size + this->next_free)
		return false;
	if (section_type == 'r' or section_type == 'm') {
		Block_section_reader * section = Block_section_reader::construct_section(this);
		section->jump_section();
		delete section;
		return true;
	}
	return false;
}


// ----- Sections -----

char Kff_file::read_section_type() {
	// Verify that header has been read.
	if (not this->header_over) {
		this->complete_header();
	}

	return (char)this->peek_byte();
}


Section::Section(Kff_file * file) {
	this->file = file;

	if (not file->header_over and file->footer_discovery_ended) {
		file->complete_header();
	}

	this->beginning = file->tellp();
}

void Section::close() {
	this->file = nullptr;
}


Section * SectionBuilder::build(Kff_file * file) {
	char type = file->read_section_type();
	switch (type) {
		case 'i':
			return new Section_Index(file);
		case 'v':
			return new Section_GV(file);
		case 'r':
			return new Section_Raw(file);
		case 'm':
			return new Section_Minimizer(file);
		default:
			cerr << "Unknown section " << type << "(" << (uint)type << ")" << endl;
			throw std::runtime_error("Unknown section type " + std::string(1, type));
	}
}


// ----- Global variables sections -----

Section_GV::Section_GV(Kff_file * file) : Section(file) {
	this->nb_vars = 0;
	this->file->global_vars.clear();

	if (this->file->is_reader) {
		this->read_section();
	}

	if (file->is_writer) {
		if (file->indexed)
			file->register_position('v');
		char type = 'v';
		this->file->write((uint8_t *)&type, 1);
	}
}

void Section_GV::write_var(const string & var_name, uint64_t value) {
	this->vars[var_name] = value;
	// A variable written twice is stored once
	this->nb_vars = this->vars.size();
	this->file->global_vars[var_name] = value;
}

void Section_GV::read_section() {
	char type = '\0';
	this->file->read((uint8_t *)&type, 1);
	if (type != 'v')
		throw "The section do not start with the 'v' char, you can't open a Global Variable section.";

	uint8_t buff[8];
	this->file->read(buff, 8);
	load_big_endian(buff, 8, this->nb_vars);
	for (uint64_t i=0 ; i<nb_vars ; i++) {
		this->read_var();
	}
}

void Section_GV::read_var() {
	if (file->tellp() >= file->end_position)
		throw "eof reached before the end of the variable section";

	// Name reading
	stringstream ss;
	char c = 'o';
	this->file->read((uint8_t *)&c, 1);
	while (c != '\0') {
		ss << c;
		this->file->read((uint8_t *)&c, 1);
	}

	// Value reading
	uint64_t value = 0;
	uint8_t buff[8];
	this->file->read(buff, 8);
	load_big_endian(buff, 8, value);

	// Saving
	string name = ss.str();
	this->vars[name] = value;
	this->file->global_vars[name] = value;
}

void Section_GV::copy(Kff_file * file) {
	// Remove empty variable sections
	if (this->vars.size() == 0)
		return;

	// Open the copy
	Section_GV sgv(file);
	// Copy all the variables
	for (const auto & it : this->vars) {
		sgv.write_var(it.first, it.second);
	}
	// Clos the copy
	sgv.close();
}

void Section_GV::close() {
	if (file->is_writer) {
		uint8_t buff[8];
		// write the number of block values
		store_big_endian(buff, 8, this->nb_vars);
		this->file->write(buff, 8);
		// Write the variables
		for (std::map<std::string,uint64_t>::iterator var_tuple=this->vars.begin() ; var_tuple != this->vars.end() ; var_tuple++) {
			const string & name = var_tuple->first;
			this->file->write((const uint8_t *)name.c_str(), name.length()+1);
			store_big_endian(buff, 8, var_tuple->second);
			this->file->write(buff, 8);
		}
	}

	Section::close();
}



Section_Index::Section_Index(Kff_file * file) : Section(file) {
	char type;
	uint8_t buff[8];

	this->next_index = 0;

	if (this->file->is_reader) {
		this->file->read((uint8_t *)&type, 1);
		if (type != 'i')
			throw "The section do not start with the 'i' char, you can not open an Index section.";

		uint64_t nb_vars;
		this->file->read(buff, 8);
		load_big_endian(buff, 8, nb_vars);
		for (uint64_t i=0 ; i<nb_vars ; i++) {
			int64_t idx = 0;
			this->file->read((uint8_t *)&type, 1);
			this->file->read(buff, 8);
			load_big_endian(buff, 8, idx);
			this->index[idx] = type;
		}

		if (nb_vars != this->index.size())
			throw "index collision in i section";

		this->file->read(buff, 8);
		load_big_endian(buff, 8, this->next_index);
	}
}

void Section_Index::register_section(char section_type, int64_t pos) {
	this->index[pos] = section_type;
}

void Section_Index::set_next_index(int64_t index) {
	this->next_index = index;
}

// void Section_Index::copy(Kff_file * file) {
// 	cerr << "You are trying to copy an index from a file to another. ";
// 	cerr << "As the positions can be different between the files, this operation is not allowed." << endl;
// 	exit(2);
// }

void Section_Index::close() {
	if (this->file->is_writer) {
		uint8_t buff[8];
		// Section header
		char type = 'i';
		this->file->write((uint8_t *)&type, 1);
		store_big_endian(buff, 8, this->index.size());
		this->file->write(buff, 8);
		// Write index
		for (std::map<int64_t, char>::iterator it=this->index.begin(); it!=this->index.end(); ++it) {
			// Section type
			type = it->second;
			this->file->write((uint8_t *)&type, 1);
			// Section index
			store_big_endian(buff, 8, it->first);
			this->file->write(buff, 8);
		}
		store_big_endian(buff, 8, this->next_index);
		this->file->write(buff, 8);
	}

	Section::close();
}



Block_section_reader * Block_section_reader::construct_section(Kff_file * file) {
	// Very and complete if needed the header
	file->complete_header();

	char type = file->read_section_type();
	if (type == 'r') {
		return new Section_Raw(file);
	} else if (type == 'm') {
		return new Section_Minimizer(file);
	} else
		return nullptr;
}


// ----- Raw sequence section -----

Section_Raw::Section_Raw(Kff_file * file) : Section(file){
	if (file->global_vars.find("k") == file->global_vars.end())
		throw "Impossible to read the raw section due to missing k variable";
	if(file->global_vars.find("max") == file->global_vars.end())
		throw "Impossible to read the raw section due to missing max variable";
	if(file->global_vars.find("data_size") == file->global_vars.end())
		throw "Impossible to read the raw section due to missing data_size variable";
	
	uint64_t k = file->global_vars["k"];
	uint64_t max = file->global_vars["max"];
	uint64_t data_size = file->global_vars["data_size"];

	this->nb_blocks = 0;

	this->k = k;
	this->max = max;
	this->data_size = data_size;

	// Computes the number of bytes needed to store the number of kmers in each block
	uint64_t nb_bits = static_cast<uint64_t>(ceil(log2(max)));
	this->nb_kmers_bytes = static_cast<uint8_t>(bytes_from_bit_array(nb_bits, 1));

	if (file->is_reader) {
		this->read_section_header();
	}

	if (file->is_writer) {
		if (file->indexed)
			file->register_position('r');
		char type = 'r';
		this->file->write((uint8_t *)&type, 1);
		this->file->write((uint8_t *)&this->nb_blocks, 8);
	}
}

uint32_t Section_Raw::read_section_header() {
	char type;
	this->file->read((uint8_t *)&type, 1);
	if (type != 'r')
		throw "The section do not start with the 'r' char, you can't open a Raw sequence section.";

	uint8_t buff[8];
	this->file->read(buff, 8);
	load_big_endian(buff, 8, this->nb_blocks);
	this->remaining_blocks = this->nb_blocks;

	return this->nb_blocks;
}

void Section_Raw::write_compacted_sequence(uint8_t* seq, uint64_t seq_size, uint8_t * data_array) {
	uint8_t buff[8];
	// 1 - Write nb kmers
	uint64_t nb_kmers = seq_size - k + 1;
	store_big_endian(buff, this->nb_kmers_bytes, nb_kmers);
	this->file->write(buff, this->nb_kmers_bytes);
	// 2 - Write sequence
	uint64_t seq_bytes_needed = (seq_size + 3) / 4;
	this->file->write(seq, seq_bytes_needed);
	// 3 - Write data
	uint64_t data_bytes_needed = data_size * nb_kmers;
	this->file->write(data_array, data_bytes_needed);

	this->nb_blocks += 1;
}

uint64_t Section_Raw::read_compacted_sequence(uint8_t* seq, uint8_t* data) {
	uint8_t buff[8];
	uint64_t nb_kmers_in_block = 1;
	// 1 - Read the number of kmers in the sequence
	if (nb_kmers_bytes != 0) {
		file->read(buff, this->nb_kmers_bytes);
		load_big_endian(buff, this->nb_kmers_bytes, nb_kmers_in_block);
	}
	// 2 - Read the sequence
	size_t seq_size = nb_kmers_in_block + k - 1;
	size_t seq_bytes_needed = (seq_size + 3) / 4;
	this->file->read(seq, seq_bytes_needed);
	// 3 - Read the data.
	uint64_t data_bytes_used = data_size * nb_kmers_in_block;
	this->file->read(data, data_bytes_used);

	this->remaining_blocks -= 1;

	return nb_kmers_in_block;
}


uint64_t Section_Raw::read_compacted_sequence(uint8_t* seq_data) {
	uint8_t buff[8];
	uint64_t nb_kmers_in_block = 1;
	// 1 - Read the number of kmers in the sequence
	if (nb_kmers_bytes != 0) {
		file->read(buff, this->nb_kmers_bytes);
		load_big_endian(buff, this->nb_kmers_bytes, nb_kmers_in_block);
	}
	// 2 - Read the sequence
	size_t seq_size = nb_kmers_in_block + k - 1;
	size_t seq_bytes_needed = (seq_size + 3) / 4;
	uint64_t data_bytes_used = data_size * nb_kmers_in_block;
	this->file->read(seq_data, seq_bytes_needed + data_bytes_used);

	this->remaining_blocks -= 1;

	return nb_kmers_in_block;
}


void Section_Raw::copy(Kff_file * file) {
	uint max_nucl = this->k + this->max - 1;
	uint8_t * seq_buffer = new uint8_t[(max_nucl + 3) / 4];
	uint8_t * data_buffer = new uint8_t[this->max * this->data_size];

	// Open the copy
	Section_Raw sr(file);
	// Copy all the variables
	for (uint i=0 ; i<this->nb_blocks ; i++) {
		// Read
		uint64_t size = this->read_compacted_sequence(seq_buffer, data_buffer);
		// Rewrite
		sr.write_compacted_sequence(seq_buffer, this->k + size - 1, data_buffer);
	}
	// Clos the copy
	sr.close();
	delete[] seq_buffer;
	delete[] data_buffer;
}


void Section_Raw::jump_sequence() {
	uint8_t buff[8];
	uint64_t nb_kmers_in_block = 1;
	// 1 - Read the number of kmers in the sequence
	if (nb_kmers_bytes != 0) {
		file->read(buff, this->nb_kmers_bytes);
		load_big_endian(buff, this->nb_kmers_bytes, nb_kmers_in_block);
	}
	// 2 - Determine the sequence size
	size_t seq_size = nb_kmers_in_block + k - 1;
	size_t seq_bytes_needed = (seq_size + 3) / 4;
	// 3 - Determine the data size
	size_t data_bytes_used = data_size * nb_kmers_in_block;
	// 4 - Jumb over the 
	file->jump(seq_bytes_needed + data_bytes_used);
	this->remaining_blocks -= 1;
}


void Section_Raw::close() {
	if (this->file->is_writer) {
		uint8_t buff[8];
		store_big_endian(buff, 8, this->nb_blocks);
		this->file->write_at(buff, 8, this->beginning + 1);
	}

	if (file->is_reader) {
		// Jump over remaining sequences of the section
		while (this->remaining_blocks > 0)
			this->jump_sequence();
	}

	Section::close();
}



// ----- Minimizer sequence section -----

Section_Minimizer::Section_Minimizer(Kff_file * file) : Section(file) {
	if (file->global_vars.find("k") == file->global_vars.end())
		throw "Impossible to read the minimizer section due to missing k variable";
	if (file->global_vars.find("m") == file->global_vars.end())
		throw "Impossible to read the minimizer section due to missing m variable";
	if(file->global_vars.find("max") == file->global_vars.end())
		throw "Impossible to read the minimizer section due to missing max variable";
	if(file->global_vars.find("data_size") == file->global_vars.end())
		throw "Impossible to read the minimizer section due to missing data_size variable";
	
	uint64_t k = file->global_vars["k"];
	uint64_t m = file->global_vars["m"];
	uint64_t max = file->global_vars["max"];
	uint64_t data_size = file->global_vars["data_size"];

	this->nb_blocks = 0;
	this->remaining_blocks = 0;

	this->k = k;
	this->m = m;
	this->max = max;
	this->data_size = data_size;

	// Computes the number of bytes needed to store the number of kmers in each block
	uint64_t nb_bits = static_cast<uint64_t>(ceil(log2(max)));
	this->nb_kmers_bytes = static_cast<uint8_t>(bytes_from_bit_array(nb_bits, 1));
	this->nb_bytes_mini = static_cast<uint8_t>(bytes_from_bit_array(2, m));
	this->minimizer = new uint8_t[nb_bytes_mini];
	memset(this->minimizer, 0, nb_bytes_mini);
	uint64_t mini_pos_bits = static_cast<uint8_t>(ceil(log2(k+max-1)));
	this->mini_pos_bytes = bytes_from_bit_array(mini_pos_bits, 1);

	if (file->is_reader) {
		this->read_section_header();
	}

	if (file->is_writer) {
		if (file->indexed)
			file->register_position('m');

		char type = 'm';
		this->file->write((uint8_t *)&type, 1);
		this->file->write(this->minimizer, this->nb_bytes_mini);
		this->file->write((uint8_t *)&this->nb_blocks, 8);
	}
}

Section_Minimizer::~Section_Minimizer() {
	delete[] this->minimizer;
}

Section_Minimizer& Section_Minimizer::operator= ( Section_Minimizer && sm) {
	file = sm.file;
	sm.file = nullptr;
	beginning = sm.beginning;
	nb_blocks = sm.nb_blocks;

	m = sm.m;
	k = sm.k;
	max = sm.max;
	data_size = sm.data_size;

	this->remaining_blocks = sm.remaining_blocks;
	this->nb_kmers_bytes = sm.nb_kmers_bytes;
	this->mini_pos_bytes = sm.mini_pos_bytes;

	nb_bytes_mini = sm.nb_bytes_mini;
	std::swap(minimizer, sm.minimizer);

	return *this;
}

void Section_Minimizer::write_minimizer(uint8_t * minimizer) {
	this->file->write_at(minimizer, this->nb_bytes_mini, this->beginning+1);

	// uint64_t pos = file->fs.tellp();
	// file->fs.seekp(this->beginning+1);
	// file->fs.write((char *)minimizer, this->nb_bytes_mini);
	// memcpy(this->minimizer, minimizer, this->nb_bytes_mini);
	// file->fs.seekp(pos);
}

uint32_t Section_Minimizer::read_section_header() {
	// Verify section type
	char type;
	this->file->read((uint8_t *)&type, 1);
	if (type != 'm')
		throw "The section do not start with the 'm' char, you can't open a Minimizer sequence section.";

	// Read the minimizer
	this->file->read(this->minimizer, this->nb_bytes_mini);

	// Read the number of following blocks
	uint8_t buff[8];
	this->file->read(buff, 8);
	load_big_endian(buff, 8, this->nb_blocks);
	this->remaining_blocks = this->nb_blocks;
	return nb_blocks;
}

void Section_Minimizer::write_compacted_sequence_without_mini(uint8_t* seq, uint64_t seq_size, uint64_t mini_pos, uint8_t * data_array) {
	uint8_t buff[8];
	// 1 - Write nb kmers
	uint64_t nb_kmers = seq_size + m - k + 1;
	store_big_endian(buff, this->nb_kmers_bytes, nb_kmers);
	file->write(buff, this->nb_kmers_bytes);
	// 2 - Write minimizer position
	store_big_endian(buff, this->mini_pos_bytes, mini_pos);
	file->write(buff, this->mini_pos_bytes);
	// 3 - Write sequence with chopped minimizer
	uint64_t seq_bytes_needed = bytes_from_bit_array(2, seq_size);
	this->file->write(seq, seq_bytes_needed);
	// 4 - Write data
	uint64_t data_bytes_needed = bytes_from_bit_array(data_size*8, nb_kmers);
	this->file->write(data_array, data_bytes_needed);

	this->nb_blocks += 1;
}

uint64_t Section_Minimizer::read_compacted_sequence_without_mini(uint8_t* seq, uint8_t* data, uint64_t & mini_pos) {
	uint8_t buff[8];
	uint64_t nb_kmers_in_block = 1;
	// 1 - Read the number of kmers in the sequence
	if (nb_kmers_bytes != 0) {
		file->read(buff, this->nb_kmers_bytes);
		load_big_endian(buff, this->nb_kmers_bytes, nb_kmers_in_block);
	}
	// 2 - Read the minimizer position
	file->read(buff, this->mini_pos_bytes);
	load_big_endian(buff, this->mini_pos_bytes, mini_pos);
	// 3 - Read the sequence
	size_t seq_size = nb_kmers_in_block + k - m - 1;
	size_t seq_bytes_needed = bytes_from_bit_array(2, seq_size);
	file->read(seq, seq_bytes_needed);
	// 4 - Read the data
	uint64_t data_bytes_needed = bytes_from_bit_array(data_size*8, nb_kmers_in_block);
	file->read(data, data_bytes_needed);
	// cout << data_bytes_needed << endl;

	this->remaining_blocks -= 1;
	return nb_kmers_in_block;
}

uint64_t Section_Minimizer::read_compacted_sequence_without_mini(uint8_t* seq_data, uint64_t & mini_pos) {
	uint8_t buff[8];
	uint64_t nb_kmers_in_block = 1;
	// 1 - Read the number of kmers in the sequence
	if (nb_kmers_bytes != 0) {
		file->read(buff, this->nb_kmers_bytes);
		load_big_endian(buff, this->nb_kmers_bytes, nb_kmers_in_block);
	}
	// 2 - Read the minimizer position
	file->read(buff, this->mini_pos_bytes);
	load_big_endian(buff, this->mini_pos_bytes, mini_pos);
	// 3 - Read the sequence + data
	size_t seq_size = nb_kmers_in_block + k - m - 1;
	size_t seq_bytes_needed = bytes_from_bit_array(2, seq_size);
	uint64_t data_bytes_needed = bytes_from_bit_array(data_size*8, nb_kmers_in_block);
	file->read(seq_data, seq_bytes_needed + data_bytes_needed);

	this->remaining_blocks -= 1;
	return nb_kmers_in_block;
}

void Section_Minimizer::copy(Kff_file * file) {
	uint max_nucl = this->k + this->max - 1;
	uint8_t * seq_buffer = new uint8_t[(max_nucl + 3) / 4];
	uint8_t * data_buffer = new uint8_t[this->max * this->data_size];
	uint64_t mini_pos = 0;

	// Open the copy
	Section_Minimizer sm(file);
	sm.write_minimizer(this->minimizer);

	// Copy all the variables
	for (uint i=0 ; i<this->nb_blocks ; i++) {
		// Read
		uint64_t size = this->read_compacted_sequence_without_mini(seq_buffer, data_buffer, mini_pos);
		// Rewrite
		sm.write_compacted_sequence_without_mini(seq_buffer, this->k + size - 1 - this->m, mini_pos, data_buffer);
	}
	// Clos the copy
	sm.close();
	delete[] seq_buffer;
	delete[] data_buffer;
}

void Section_Minimizer::jump_sequence() {
	uint64_t nb_kmers_in_block = 1;
	uint8_t buff[8];
	// 1 - Read the number of kmers in the sequence
	if (nb_kmers_bytes != 0) {
		file->read(buff, this->nb_kmers_bytes);
		// 2 - Convert from big endian
		load_big_endian(buff, this->nb_kmers_bytes, nb_kmers_in_block);
	}
	// 3 - Determine the sequence size
	size_t seq_size = nb_kmers_in_block + k - m - 1;
	size_t seq_bytes_needed = (seq_size + 3) / 4;
	// 3 - Determine the data size
	size_t data_bytes_used = data_size * nb_kmers_in_block;
	// 4 - Jumb over the 
	file->jump((long)(this->mini_pos_bytes + seq_bytes_needed + data_bytes_used));

	this->remaining_blocks -= 1;
}



/* Similar to the previous function but on the right */
static void rightshift8(uint8_t * bitarray, size_t length, size_t bitshift) {
	assert(bitshift < 8);

	if (length > 0) {
		for (uint64_t i=length-1 ; i>0 ; i--) {
			bitarray[i] = (bitarray[i-1] << (8-bitshift)) | (bitarray[i] >> bitshift);
		}
		bitarray[0] >>= bitshift;
	}
}


void Section_Minimizer::write_compacted_sequence (uint8_t* seq, uint64_t seq_size, uint64_t mini_pos, uint8_t * data_array) {
	// The sequence is copied without its minimizer, right aligned with null padding bits
	const uint64_t no_mini_size = seq_size - m;
	const uint64_t seq_offset = (4 - seq_size % 4) % 4;
	const uint64_t no_mini_offset = (4 - no_mini_size % 4) % 4;
	this->nucleotide_buffer.assign(bytes_from_bit_array(2, no_mini_size) + 1, 0);
	uint8_t * no_mini = this->nucleotide_buffer.data();

	// Prefix then suffix of the minimizer
	copy_nucleotides(no_mini, no_mini_offset, seq, seq_offset, mini_pos);
	copy_nucleotides(no_mini, no_mini_offset + mini_pos, seq, seq_offset + mini_pos + m, seq_size - mini_pos - m);

	// Write the compacted sequence into file.
	this->write_compacted_sequence_without_mini(no_mini, no_mini_size, mini_pos, data_array);
}


void Section_Minimizer::add_minimizer(uint64_t nb_kmer, uint8_t * seq, uint64_t mini_pos) {
	const uint64_t seq_size = nb_kmer + k - 1;
	const uint64_t seq_bytes = bytes_from_bit_array(2, seq_size);
	const uint64_t no_mini_size = seq_size - m;
	// Padding nucleotides of the result, of the input and of the minimizer
	const uint64_t seq_offset = (4 - seq_size % 4) % 4;
	const uint64_t no_mini_offset = (4 - no_mini_size % 4) % 4;
	const uint64_t mini_offset = (4 - m % 4) % 4;

	// The full sequence is built right aligned with null padding bits
	this->nucleotide_buffer.assign(seq_bytes, 0);
	uint8_t * full = this->nucleotide_buffer.data();
	copy_nucleotides(full, seq_offset, seq, no_mini_offset, mini_pos);
	copy_nucleotides(full, seq_offset + mini_pos, this->minimizer, mini_offset, m);
	copy_nucleotides(full, seq_offset + mini_pos + m, seq, no_mini_offset + mini_pos, no_mini_size - mini_pos);
	memcpy(seq, full, seq_bytes);
}



uint64_t Section_Minimizer::read_compacted_sequence(uint8_t* seq, uint8_t* data) {
	// Read the block
	uint64_t mini_pos;
	uint64_t nb_kmers_in_block = this->read_compacted_sequence_without_mini(seq, data, mini_pos);
	this->add_minimizer(nb_kmers_in_block, seq, mini_pos);
	
	return nb_kmers_in_block;
}

uint64_t Section_Minimizer::read_compacted_sequence(uint8_t* seq_data) {
	// Read the block
	uint64_t mini_pos;
	uint64_t nb_kmers_in_block = this->read_compacted_sequence_without_mini(seq_data, mini_pos);
	
	// Determine the number of new bytes needed for minimizer insertion
	uint64_t seq_size_nucls = k - m + nb_kmers_in_block - 1;
	uint64_t free_nucls = (4 - seq_size_nucls) % 4;
	uint64_t bytes_needed = (m - free_nucls + 3) / 4;
	
	if (bytes_needed > 0)
	{
		// Move the data to the right to have the hole needed for minimizer insertion
		uint64_t seq_size_bytes = (seq_size_nucls + 3) / 4;
		uint64_t data_size_bytes = data_size * nb_kmers_in_block;
		for (size_t i=0 ; i<data_size_bytes ; i++)
		{
			// Shift the ith byte to "byte_needed" bytes on the right
			size_t byte_idx = seq_size_bytes + data_size_bytes - 1 - i;
			seq_data[byte_idx + bytes_needed] = seq_data[byte_idx];
			seq_data[byte_idx] = 0;
		}
	}

	// Insert minimizer
	this->add_minimizer(nb_kmers_in_block, seq_data, mini_pos);
	
	return nb_kmers_in_block;
}

void Section_Minimizer::close() {
	if (file->is_writer) {
		uint8_t tmp[8];
		store_big_endian(tmp, 8, this->nb_blocks);
		this->file->write_at(tmp, 8, this->beginning + 1l + (long)this->nb_bytes_mini);
	}

	if (file->is_reader) {
		// Jump over remaining sequences of the section
		while (this->remaining_blocks > 0)
			this->jump_sequence();
	}

	Section::close();
}





// -------- Start of the high level API -----------

Kff_reader::Kff_reader(std::string filename) {
	// Open the file
	this->file = new Kff_file(filename, "r");

	this->current_seq_data = new uint8_t[1];
	this->current_seq_data[0] = 0;

	// Create fake small datastrucutes waiting for the right values.
	this->current_shifts = new uint8_t*[4];
	this->current_shifts[0] = this->current_seq_data;
	for (uint8_t i=1 ; i<4 ; i++) {
		this->current_shifts[i] = new uint8_t[1];
		this->current_shifts[i][0] = 0;
	}

	this->current_section = NULL;
	this->current_kmer = new uint8_t[1];
	this->remaining_kmers = 0;

	this->k = 0;
	this->max = 0;
	this->data_size = 0;

	this->has_next();
}

Kff_reader::~Kff_reader() {
	delete[] this->current_kmer;
	delete[] this->current_seq_data;
	if (this->current_section != NULL)
		delete this->current_section;

	for (uint i=1 ; i<4 ; i++)
		delete[] this->current_shifts[i];
	delete[] this->current_shifts;

	delete this->file;
}

void Kff_reader::read_until_first_section_block() {
	while (current_section == NULL or remaining_blocks == 0) {
		if (this->file->tellp() == this->file->end_position) {
			break;
		}

		// char section_type = this->file->read_section_type();
		char section_type = file->read_section_type();
		// --- Update data structure sizes ---
		if (section_type == 'v')
		{
			// Read the global variable block
			Section_GV gvs(file);
			// Update sequence size if k or max change
			if (gvs.vars.find("k") != gvs.vars.end()
				or gvs.vars.find("max") != gvs.vars.end())
			{
				// Compute the max size of a sequence
				this->k = this->file->global_vars["k"];
				this->max = this->file->global_vars["max"];
				uint64_t seq_max_size = bytes_from_bit_array(2, max + k - 1);
				uint64_t data_max_size = data_size * max;
				// sequence + data buffer
				delete[] this->current_seq_data;
				this->current_seq_data = new uint8_t[seq_max_size + data_max_size];
				
				// Shifts
				this->current_shifts[0] = this->current_seq_data;
				for (uint8_t i=1 ; i<4 ; i++)
				{
					delete[] this->current_shifts[i];
					this->current_shifts[i] = new uint8_t[seq_max_size];
					memset(this->current_shifts[i], 0, seq_max_size);
				}

				// Current kmer
				delete[] this->current_kmer;
				this->current_kmer = new uint8_t[k/4 + 1];
				memset(this->current_kmer, 0, (k/4+1));
			}

			// Update the data array size
			if (gvs.vars.find("data_size") != gvs.vars.end()
				or gvs.vars.find("max") != gvs.vars.end())
			{
				// Update global variables
				this->max = this->file->global_vars["max"];
				this->data_size = this->file->global_vars["data_size"];
				// Update buffer size
				uint64_t seq_max_size = bytes_from_bit_array(2, max + k - 1);
				uint64_t data_max_size = this->data_size * max;
				delete[] this->current_seq_data;
				this->current_seq_data = new uint8_t[seq_max_size + data_max_size];
				this->current_shifts[0] = this->current_seq_data;
				memset(this->current_seq_data, 0, seq_max_size + data_max_size);
			}
		}
		// Mount data from the files to the datastructures.
		else if (section_type == 'i') {
			Section_Index index(file);
			index.close();
		} else {
			current_section = Block_section_reader::construct_section(file);
			remaining_blocks = current_section->nb_blocks;
			// Sections without block are skipped
			if (remaining_blocks == 0) {
				delete current_section;
				current_section = NULL;
			}
		}
	}
}


void Kff_reader::read_next_block() {
	// Read from the file
	current_seq_kmers = remaining_kmers = current_section->read_compacted_sequence(current_seq_data);
	current_seq_nucleotides = remaining_kmers + current_section->k - 1;
	current_seq_bytes = bytes_from_bit_array(2, current_seq_nucleotides);

	// Create the 4 possible shifts of the sequence for easy use.
	for (uint8_t i=1 ; i<min((uint64_t)4, remaining_kmers) ; i++) {
		// Copy
		memcpy(current_shifts[i], current_shifts[0], current_seq_bytes);
		// Shift
		rightshift8(current_shifts[i], current_seq_bytes, 2 * i);
	}
}

bool Kff_reader::has_next() {
	if (current_section == NULL and (file->end_position > file->tellp()))
		read_until_first_section_block();
	return file->end_position > file->tellp();
}

uint64_t Kff_reader::next_block(uint8_t* & sequence, uint8_t* & data) {
	// Verify the abylity to find another kmer in the file.
	if (!this->has_next()){
		sequence = NULL;
		data = NULL;
		return 0;
	}

	uint64_t nb_kmers = current_section->read_compacted_sequence(sequence, data);
	
	remaining_kmers = 0;
	remaining_blocks -= 1;
	if (remaining_blocks == 0) {
		delete current_section;
		current_section = NULL;
	}

	return nb_kmers;
}

bool Kff_reader::next_kmer(uint8_t* & kmer, uint8_t* & data) {
	// Verify the abylity to find another kmer in the file.
	if (!this->has_next()){
		kmer = NULL;
		data = NULL;
		return false;
	}

	// Load the next block
	if (remaining_kmers == 0) {
		read_next_block();
	}

	uint64_t right_shift = (remaining_kmers - 1) % 4;
	uint64_t prefix_offset = (4 - (current_seq_nucleotides % 4)) % 4;
	uint64_t kmer_idx = current_seq_kmers - remaining_kmers;

	uint64_t start_nucl = prefix_offset + right_shift + kmer_idx;
	uint64_t start_byte = start_nucl / 4;
	uint64_t end_nucl = start_nucl + this->k - 1;
	uint64_t end_byte = end_nucl / 4;

	memcpy(current_kmer, current_shifts[right_shift]+start_byte, end_byte-start_byte+1);
	kmer = current_kmer;
	data = current_seq_data + current_seq_bytes + (current_seq_kmers - remaining_kmers) * this->data_size;
	
	// Read the next block if needed.
	remaining_kmers -= 1;
	if (remaining_kmers == 0) {
		remaining_blocks -= 1;
		if (remaining_blocks == 0) {
			delete current_section;
			current_section = NULL;
		}
	}

	return true;
}


uint64_t Kff_reader::get_var(string name) {
	if (file->global_vars.find(name) != file->global_vars.end())
		return file->global_vars[name];

	cerr << "Variable " << name << " is absent from the file." << endl;
	exit(2);

	return 0;
}


uint8_t * Kff_reader::get_encoding() {
	return file->encoding;
}












