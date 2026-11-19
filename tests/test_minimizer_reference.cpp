/**
 * Differential tests: the minimizer insertion (reading) and removal (writing)
 * of the library are compared to a verbatim copy of the original
 * implementation (kff-cpp-api 49dd236) on random inputs, including random
 * bits in the padding of the inputs and canary bytes after the buffers.
 * The library must produce exactly the same bytes.
 */
#include <cstring>

#include "framework.hpp"
#include "helpers.hpp"

using namespace kfftest;

namespace reference {

static uint64_t bytes_from_bit_array(uint64_t bits_per_elem, uint64_t nb_elem) {
	if (bits_per_elem == 0 or nb_elem == 0)
		return 0;
	else
		return ((bits_per_elem * nb_elem - 1) / 8) + 1;
}

static void leftshift8(uint8_t * bitarray, size_t length, size_t bitshift) {
	if (length > 0) {
		for (uint64_t i=0 ; i<length-1 ; i++) {
			bitarray[i] = (bitarray[i] << bitshift) | (bitarray[i+1] >> (8-bitshift));
		}
		bitarray[length-1] <<= bitshift;
	}
}

static void rightshift8(uint8_t * bitarray, size_t length, size_t bitshift) {
	if (length > 0) {
		for (uint64_t i=length-1 ; i>0 ; i--) {
			bitarray[i] = (bitarray[i-1] << (8-bitshift)) | (bitarray[i] >> bitshift);
		}
		bitarray[0] >>= bitshift;
	}
}

static uint8_t fusion8(uint8_t left_bits, uint8_t right_bits, size_t merge_index) {
	uint8_t mask = 0xFF << (8-merge_index);
	return (left_bits & mask) | (right_bits & ~mask);
}

/** Section_Minimizer::add_minimizer */
static void add_minimizer(uint64_t k, uint64_t m, const uint8_t * minimizer, uint64_t nb_kmer, uint8_t * seq, uint64_t mini_pos) {
	const uint64_t nb_bytes_mini = bytes_from_bit_array(2, m);
	uint64_t seq_size = nb_kmer + k - 1;
	uint64_t seq_bytes = bytes_from_bit_array(2, seq_size);
	uint64_t seq_left_offset = (4 - (seq_size % 4)) % 4;
	uint64_t no_mini_size = seq_size - m;
	uint64_t no_mini_bytes = bytes_from_bit_array(2, no_mini_size);
	uint64_t no_mini_left_offset = (4 - (no_mini_size % 4)) % 4;
	leftshift8(seq, no_mini_bytes, no_mini_left_offset*2);

	uint8_t * suffix = new uint8_t[seq_bytes];
	memset(suffix, 0, seq_bytes);
	unsigned suff_nucl = seq_size - m - mini_pos;
	unsigned no_mini_suff_start_nucl = mini_pos;
	unsigned no_mini_suff_start_byte = no_mini_suff_start_nucl / 4;
	unsigned no_mini_suff_bytes = no_mini_bytes - no_mini_suff_start_byte;
	memcpy(suffix, seq + no_mini_suff_start_byte, no_mini_suff_bytes);
	unsigned no_mini_suff_offset = no_mini_suff_start_nucl % 4;
	leftshift8(suffix, no_mini_suff_bytes, no_mini_suff_offset * 2);

	uint8_t * mini = new uint8_t[seq_bytes];
	memset(mini, 0, seq_bytes);
	memcpy(mini, minimizer, nb_bytes_mini);
	unsigned mini_offset = (4 - (m % 4)) % 4;
	leftshift8(mini, nb_bytes_mini, mini_offset * 2);

	unsigned final_mini_start_nucl = mini_pos;
	unsigned final_mini_start_byte = mini_pos / 4;
	unsigned final_mini_offset = final_mini_start_nucl % 4;
	unsigned final_mini_byte_size = (m + final_mini_offset + 3) / 4;
	rightshift8(mini, seq_bytes, final_mini_offset * 2);

	seq[final_mini_start_byte] = fusion8(seq[final_mini_start_byte], mini[0], final_mini_offset * 2);
	for (unsigned idx=1 ; idx<final_mini_byte_size ; idx++) {
		seq[final_mini_start_byte+idx] = mini[idx];
	}

	unsigned final_suff_start_nucl = final_mini_start_nucl + m;
	unsigned final_suff_start_byte = final_suff_start_nucl / 4;
	unsigned final_suff_offset = final_suff_start_nucl % 4;
	unsigned final_suff_byte_size = (suff_nucl + final_suff_offset + 3) / 4;
	if (final_suff_byte_size > 0)
	{
		rightshift8(suffix, seq_bytes, final_suff_offset * 2);
		seq[final_suff_start_byte] = fusion8(seq[final_suff_start_byte], suffix[0], final_suff_offset * 2);
		for (uint64_t idx=1 ; idx<final_suff_byte_size ; idx++) {
			seq[final_suff_start_byte+idx] = suffix[idx];
		}
	}

	rightshift8(seq, seq_bytes, seq_left_offset * 2);

	delete[] suffix;
	delete[] mini;
}

/** Sequence transformation of Section_Minimizer::write_compacted_sequence.
 *  Returns the bytes given to write_compacted_sequence_without_mini. */
static Bytes remove_minimizer(uint64_t m, const uint8_t * seq, uint64_t seq_size, uint64_t mini_pos) {
	uint64_t seq_bytes = bytes_from_bit_array(2, seq_size);
	unsigned left_offset_nucl = (4 - seq_size % 4) % 4;

	uint8_t * seq_copy = new uint8_t[seq_bytes];
	memcpy(seq_copy, seq, seq_bytes);

	unsigned mini_start_byte = (mini_pos + left_offset_nucl) / 4;
	unsigned suff_start_byte = (mini_pos + m + left_offset_nucl) / 4;
	unsigned suff_bytes = seq_bytes-suff_start_byte;
	for (unsigned i=0 ; i<suff_bytes ; i++) {
		seq_copy[mini_start_byte + i] = seq_copy[suff_start_byte + i];
	}

	unsigned mini_offset = (mini_pos + left_offset_nucl) % 4;
	unsigned suff_offset = (mini_pos + m + left_offset_nucl) % 4;
	if (mini_offset < suff_offset)
		leftshift8(seq_copy + mini_start_byte, seq_bytes - mini_start_byte, (suff_offset - mini_offset) * 2);
	else
		rightshift8(seq_copy + mini_start_byte, seq_bytes - mini_start_byte, (mini_offset - suff_offset) * 2);

	seq_copy[mini_start_byte] = fusion8(seq[mini_start_byte], seq_copy[mini_start_byte], mini_offset * 2);

	leftshift8(seq_copy, seq_bytes, left_offset_nucl * 2);
	rightshift8(seq_copy, seq_bytes, ((4 - ((seq_size - m) % 4)) % 4) * 2);

	Bytes result(seq_copy, seq_copy + bytes_from_bit_array(2, seq_size - m));
	delete[] seq_copy;
	return result;
}

} // namespace reference

/** Random bits in the padding of the first byte of a sequence of size nucleotides */
static void dirty_padding(Bytes & bytes, size_t size, Random & rnd) {
	size_t padding = (4 - size % 4) % 4;
	if (padding > 0 and not bytes.empty())
		bytes[0] |= rnd.between(0, 255) & ~((1u << (2 * (4 - padding))) - 1) & 0xff;
}

TEST_CASE("minimizer reference: add_minimizer gives the same bytes as the original code") {
	TempDir dir;
	Random rnd(100);
	Kff_file file(dir.file("scratch.kff"), "w");
	for (uint64_t m = 1; m <= 31; m++) {
		for (int config = 0; config < 6; config++) {
			const uint64_t k = m + rnd.between(0, config < 3 ? 6 : 40);
			const uint64_t max = rnd.between(1, 70);
			Section_GV sgv(&file);
			sgv.write_var("k", k);
			sgv.write_var("m", m);
			sgv.write_var("max", max);
			sgv.write_var("data_size", 0);
			sgv.close();
			Section_Minimizer sm(&file);

			for (int test = 0; test < 40; test++) {
				// Minimizer with random padding bits
				Bytes minimizer = encode(rnd.sequence(m));
				dirty_padding(minimizer, m, rnd);
				memcpy(sm.minimizer, minimizer.data(), minimizer.size());

				const uint64_t nb_kmers = rnd.between(1, max);
				const uint64_t seq_size = nb_kmers + k - 1;
				const uint64_t mini_pos = rnd.between(0, seq_size - m);
				// Sequence without minimizer, random padding and random bytes after it
				Bytes buffer = rnd.bytes((seq_size + 3) / 4 + 8);
				Bytes without = encode(rnd.sequence(seq_size - m));
				dirty_padding(without, seq_size - m, rnd);
				std::copy(without.begin(), without.end(), buffer.begin());

				Bytes expected = buffer;
				reference::add_minimizer(k, m, minimizer.data(), nb_kmers, expected.data(), mini_pos);
				sm.add_minimizer(nb_kmers, buffer.data(), mini_pos);
				CHECK_EQ(buffer, expected);
			}
			sm.close();
		}
	}
	file.close(false);
}

TEST_CASE("minimizer reference: write_compacted_sequence writes the same bytes as the original code") {
	TempDir dir;
	Random rnd(101);
	for (uint64_t m = 1; m <= 31; m++) {
		const uint64_t k = m + rnd.between(0, m % 2 ? 5 : 40);
		const uint64_t max = rnd.between(1, 70);
		const uint64_t data_size = rnd.between(0, 2);
		std::string name = dir.file("remove.kff");
		Bytes expected_section = {'m'};
		long beginning = 0;
		{
			Kff_file file(name, "w");
			file.set_indexation(false);
			Section_GV sgv(&file);
			sgv.write_var("k", k);
			sgv.write_var("m", m);
			sgv.write_var("max", max);
			sgv.write_var("data_size", data_size);
			sgv.close();

			Section_Minimizer sm(&file);
			beginning = sm.beginning;
			Bytes minimizer = encode(rnd.sequence(m));
			sm.write_minimizer(minimizer.data());
			expected_section += minimizer;
			const int nb_blocks = 50;
			expected_section += big_endian(nb_blocks, 8);
			for (int b = 0; b < nb_blocks; b++) {
				const uint64_t nb_kmers = rnd.between(1, max);
				const uint64_t seq_size = nb_kmers + k - 1;
				const uint64_t mini_pos = rnd.between(0, seq_size - m);
				Bytes seq = encode(rnd.sequence(seq_size));
				dirty_padding(seq, seq_size, rnd);
				Bytes data = rnd.bytes(nb_kmers * data_size);

				expected_section += big_endian(nb_kmers, bytes_for_max(max));
				expected_section += big_endian(mini_pos, mini_pos_bytes(k, max));
				expected_section += reference::remove_minimizer(m, seq.data(), seq_size, mini_pos);
				expected_section += data;

				Bytes seq_copy = seq;
				sm.write_compacted_sequence(seq.data(), seq_size, mini_pos, data.data());
				CHECK(seq == seq_copy);  // the input is not modified
			}
			sm.close();
			file.close();
		}
		Bytes bytes = read_file(name);
		CHECK_EQ(slice(bytes, beginning, bytes.size() - 3), expected_section);
	}
}
