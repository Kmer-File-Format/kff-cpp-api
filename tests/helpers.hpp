/**
 * Helpers for the unit tests: temporary files, an independent 2-bit
 * encoder/decoder and serializers that build the expected bytes directly
 * from the format specification (kff-reference/README.md).
 */
#ifndef KFF_TEST_HELPERS_HPP
#define KFF_TEST_HELPERS_HPP

#include <atomic>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include <unistd.h>

#include "kff_io.hpp"

namespace kfftest {

namespace fs = std::filesystem;
using Bytes = std::vector<uint8_t>;

/** A directory deleted at the end of the test */
class TempDir {
public:
	fs::path path;

	TempDir() {
		static std::atomic<int> counter{0};
		path = fs::temp_directory_path() /
		       ("kff-tests-" + std::to_string(getpid()) + "-" + std::to_string(counter++));
		fs::remove_all(path);
		fs::create_directories(path);
	}
	~TempDir() {
		std::error_code ec;
		fs::remove_all(path, ec);
	}
	std::string file(const std::string & name) const { return (path / name).string(); }
};

inline Bytes read_file(const std::string & filename) {
	std::ifstream is(filename, std::ios::binary);
	if (!is)
		throw std::runtime_error("cannot open " + filename);
	return Bytes(std::istreambuf_iterator<char>(is), std::istreambuf_iterator<char>());
}

inline void write_file(const std::string & filename, const Bytes & bytes) {
	std::ofstream os(filename, std::ios::binary);
	os.write(reinterpret_cast<const char *>(bytes.data()), bytes.size());
}

inline Bytes slice(const Bytes & bytes, size_t begin, size_t end) {
	return Bytes(bytes.begin() + begin, bytes.begin() + end);
}

/** Concatenate byte arrays */
inline Bytes & operator+=(Bytes & left, const Bytes & right) {
	left.insert(left.end(), right.begin(), right.end());
	return left;
}

/** Big endian representation of value on nb_bytes bytes */
inline Bytes big_endian(uint64_t value, size_t nb_bytes) {
	Bytes bytes(nb_bytes, 0);
	for (size_t i = 0; i < nb_bytes && i < 8; i++)
		bytes[nb_bytes - 1 - i] = static_cast<uint8_t>(value >> (8 * i));
	return bytes;
}

inline Bytes c_string(const std::string & s) {
	Bytes bytes(s.begin(), s.end());
	bytes.push_back(0);
	return bytes;
}

// ----- Nucleotides -----

/** Default encoding of the API: A=0 C=1 G=3 T=2 */
const uint8_t DEFAULT_ENCODING[4] = {0, 1, 3, 2};

inline uint8_t nucleotide_value(char nucleotide, const uint8_t encoding[4]) {
	switch (nucleotide) {
		case 'A': return encoding[0];
		case 'C': return encoding[1];
		case 'G': return encoding[2];
		case 'T': return encoding[3];
	}
	throw std::invalid_argument(std::string("not a nucleotide: ") + nucleotide);
}

/** 2 bits per nucleotide, big endian, right aligned (padding bits set to 0 on the left) */
inline Bytes encode(const std::string & sequence, const uint8_t encoding[4] = DEFAULT_ENCODING) {
	const size_t size = sequence.size();
	const size_t padding = (4 - size % 4) % 4;
	Bytes bytes((size + 3) / 4, 0);
	for (size_t i = 0; i < size; i++) {
		size_t position = padding + i;
		bytes[position / 4] |= nucleotide_value(sequence[i], encoding) << (2 * (3 - position % 4));
	}
	return bytes;
}

/** Inverse of encode. The padding bits are ignored. */
inline std::string decode(const uint8_t * bytes, size_t size, const uint8_t encoding[4] = DEFAULT_ENCODING) {
	char letters[4];
	letters[encoding[0]] = 'A';
	letters[encoding[1]] = 'C';
	letters[encoding[2]] = 'G';
	letters[encoding[3]] = 'T';
	const size_t padding = (4 - size % 4) % 4;
	std::string sequence(size, '?');
	for (size_t i = 0; i < size; i++) {
		size_t position = padding + i;
		sequence[i] = letters[(bytes[position / 4] >> (2 * (3 - position % 4))) & 0b11];
	}
	return sequence;
}

inline std::string decode(const Bytes & bytes, size_t size, const uint8_t encoding[4] = DEFAULT_ENCODING) {
	return decode(bytes.data(), size, encoding);
}

/** Deterministic pseudo-random generator for the tests */
class Random {
public:
	std::mt19937_64 engine;
	explicit Random(uint64_t seed) : engine(seed) {}

	uint64_t between(uint64_t min, uint64_t max) {  // inclusive bounds
		return std::uniform_int_distribution<uint64_t>(min, max)(engine);
	}
	std::string sequence(size_t size) {
		static const char letters[] = "ACGT";
		std::string s(size, 'A');
		for (char & c : s)
			c = letters[between(0, 3)];
		return s;
	}
	Bytes bytes(size_t size) {
		Bytes b(size);
		for (uint8_t & byte : b)
			byte = static_cast<uint8_t>(between(0, 255));
		return b;
	}
};

// ----- Expected bytes of the format elements -----

/** Number of bytes used for a value that can go up to max (ceil(log2(max)) bits) */
inline size_t bytes_for_max(uint64_t max) {
	size_t bits = 0;
	while (bits < 64 && (uint64_t(1) << bits) < max)
		bits += 1;
	return (bits + 7) / 8;
}

inline Bytes header_bytes(const uint8_t encoding[4] = DEFAULT_ENCODING, bool unique = false,
                          bool canonical = false, const Bytes & metadata = Bytes()) {
	Bytes b = {'K', 'F', 'F', KFF_VERSION_MAJOR, KFF_VERSION_MINOR,
	           static_cast<uint8_t>((encoding[0] << 6) | (encoding[1] << 4) | (encoding[2] << 2) | encoding[3]),
	           static_cast<uint8_t>(unique ? 1 : 0), static_cast<uint8_t>(canonical ? 1 : 0)};
	b += big_endian(metadata.size(), 4);
	b += metadata;
	return b;
}

/** Variables are written sorted by name */
inline Bytes variables_bytes(const std::vector<std::pair<std::string, uint64_t>> & sorted_vars) {
	Bytes b = {'v'};
	b += big_endian(sorted_vars.size(), 8);
	for (const auto & var : sorted_vars) {
		b += c_string(var.first);
		b += big_endian(var.second, 8);
	}
	return b;
}

struct Block {
	std::string sequence;  // full sequence (with the minimizer for 'm' sections)
	Bytes data;            // data_size bytes per kmer
	uint64_t mini_pos = 0; // position of the minimizer ('m' sections only)
};

inline Bytes raw_section_bytes(const std::vector<Block> & blocks, uint64_t k, uint64_t max,
                               const uint8_t encoding[4] = DEFAULT_ENCODING) {
	Bytes b = {'r'};
	b += big_endian(blocks.size(), 8);
	for (const Block & block : blocks) {
		b += big_endian(block.sequence.size() - k + 1, bytes_for_max(max));
		b += encode(block.sequence, encoding);
		b += block.data;
	}
	return b;
}

/** Width of the minimizer position field used by the implementation:
 *  ceil(log2(k + max - 1)) bits, rounded to bytes. */
inline size_t mini_pos_bytes(uint64_t k, uint64_t max) { return bytes_for_max(k + max - 1); }

inline Bytes minimizer_section_bytes(const std::string & minimizer, const std::vector<Block> & blocks,
                                     uint64_t k, uint64_t max,
                                     const uint8_t encoding[4] = DEFAULT_ENCODING) {
	const size_t m = minimizer.size();
	Bytes b = {'m'};
	b += encode(minimizer, encoding);
	b += big_endian(blocks.size(), 8);
	for (const Block & block : blocks) {
		std::string without = block.sequence.substr(0, block.mini_pos) + block.sequence.substr(block.mini_pos + m);
		b += big_endian(block.sequence.size() - k + 1, bytes_for_max(max));
		b += big_endian(block.mini_pos, mini_pos_bytes(k, max));
		b += encode(without, encoding);
		b += block.data;
	}
	return b;
}

/** Random block of nb_kmers kmers. For 'm' sections the minimizer is planted at a random position. */
inline Block random_block(Random & rnd, uint64_t k, uint64_t nb_kmers, uint64_t data_size,
                          const std::string & minimizer = "") {
	Block block;
	block.sequence = rnd.sequence(k + nb_kmers - 1);
	block.data = rnd.bytes(data_size * nb_kmers);
	if (!minimizer.empty()) {
		block.mini_pos = rnd.between(0, block.sequence.size() - minimizer.size());
		block.sequence.replace(block.mini_pos, minimizer.size(), minimizer);
	}
	return block;
}

/** Position of the first byte after the header (fixed size header + metadata) */
inline size_t first_section_position(size_t metadata_size = 0) { return 12 + metadata_size; }

} // namespace kfftest

#endif
