#include "framework.hpp"
#include "helpers.hpp"

using namespace kfftest;

/* The blocks of the sections are written in the buffer of the file when they fit in it, through
 * Kff_file::write() otherwise (the buffer is then enlarged or saved on disk). These tests write
 * many blocks of random sizes, with the usual buffer (enlarged up to 1 MB) and with a buffer that
 * is never enlarged (saved on disk each time it is full, blocks larger than the buffer), and
 * compare the bytes of the sections with the format. */

/** The write buffer of the file keeps its initial size: it is saved on disk when it is full */
static void never_enlarge_buffer(Kff_file & file) {
	file.max_buffer_size = file.buffer_size;
}

/** Random padding bits on the left of a right aligned sequence of size nucleotides */
static void dirty_padding(Bytes & bytes, size_t size, Random & rnd) {
	size_t padding = (4 - size % 4) % 4;
	if (padding > 0 and not bytes.empty())
		bytes[0] |= rnd.between(0, 255) & ~((1u << (2 * (4 - padding))) - 1) & 0xff;
}

static void write_variables(Kff_file & file, uint64_t k, uint64_t m, uint64_t max, uint64_t data_size) {
	Section_GV sgv(&file);
	sgv.write_var("k", k);
	if (m > 0)
		sgv.write_var("m", m);
	sgv.write_var("max", max);
	sgv.write_var("data_size", data_size);
	sgv.close();
}

/** Sizes of blocks: mostly small, sometimes the largest possible */
static uint64_t random_nb_kmers(Random & rnd, uint64_t max) {
	return rnd.between(0, 9) == 0 ? max : rnd.between(1, std::min<uint64_t>(max, 40));
}

TEST_CASE("block writes: raw blocks around the limits of the file buffer") {
	TempDir dir;
	Random rnd(300);
	int file_idx = 0;
	for (bool limited : {false, true})
		for (uint64_t k : {1, 5, 31, 64})
			for (uint64_t max : {1, 3, 255, 3000})
				for (uint64_t data_size : {0, 1, 3}) {
					std::string name = dir.file("raw" + std::to_string(file_idx++) + ".kff");
					std::vector<Block> blocks;
					for (int b = 0; b < 120; b++)
						blocks.push_back(random_block(rnd, k, random_nb_kmers(rnd, max), data_size));

					Kff_file file(name, "w");
					file.set_indexation(false);
					if (limited)
						never_enlarge_buffer(file);
					write_variables(file, k, 0, max, data_size);
					Section_Raw sr(&file);
					long beginning = sr.beginning;
					for (const Block & block : blocks) {
						Bytes seq = encode(block.sequence);
						Bytes data = block.data;
						sr.write_compacted_sequence(seq.data(), block.sequence.size(), data.data());
					}
					CHECK_EQ(sr.nb_blocks, (uint64_t)blocks.size());
					CHECK_EQ(file.tellp(), (unsigned long)(beginning + raw_section_bytes(blocks, k, max).size()));
					sr.close();
					file.close();

					Bytes bytes = read_file(name);
					CHECK_EQ(slice(bytes, beginning, bytes.size() - 3), raw_section_bytes(blocks, k, max));
				}
}

TEST_CASE("block writes: minimizer blocks around the limits of the file buffer") {
	TempDir dir;
	Random rnd(301);
	int file_idx = 0;
	for (bool limited : {false, true})
		for (bool without_mini : {false, true})
			for (uint64_t m : {1, 7, 13})
				for (uint64_t k : {m, m + 5, (uint64_t)31, (uint64_t)70})
					for (uint64_t max : {1, 2, 255, 2000})
						for (uint64_t data_size : {0, 2}) {
							std::string name = dir.file("mini" + std::to_string(file_idx++) + ".kff");
							std::string minimizer = rnd.sequence(m);
							std::vector<Block> blocks;
							for (int b = 0; b < 80; b++)
								blocks.push_back(random_block(rnd, k, random_nb_kmers(rnd, max), data_size, minimizer));

							Kff_file file(name, "w");
							file.set_indexation(false);
							if (limited)
								never_enlarge_buffer(file);
							write_variables(file, k, m, max, data_size);
							Section_Minimizer sm(&file);
							long beginning = sm.beginning;
							Bytes mini = encode(minimizer);
							sm.write_minimizer(mini.data());
							for (const Block & block : blocks) {
								Bytes data = block.data;
								if (without_mini) {
									std::string without = block.sequence.substr(0, block.mini_pos) + block.sequence.substr(block.mini_pos + m);
									Bytes seq = encode(without);
									sm.write_compacted_sequence_without_mini(seq.data(), without.size(), block.mini_pos, data.data());
								} else {
									Bytes seq = encode(block.sequence);
									dirty_padding(seq, block.sequence.size(), rnd);
									Bytes seq_copy = seq;
									sm.write_compacted_sequence(seq.data(), block.sequence.size(), block.mini_pos, data.data());
									CHECK(seq == seq_copy);  // the input is not modified
								}
							}
							CHECK_EQ(sm.nb_blocks, (uint64_t)blocks.size());
							Bytes expected = minimizer_section_bytes(minimizer, blocks, k, max);
							CHECK_EQ(file.tellp(), (unsigned long)(beginning + expected.size()));
							sm.close();
							file.close();

							Bytes bytes = read_file(name);
							CHECK_EQ(slice(bytes, beginning, bytes.size() - 3), expected);
						}
}
