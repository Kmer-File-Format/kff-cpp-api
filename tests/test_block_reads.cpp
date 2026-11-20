#include "framework.hpp"
#include "helpers.hpp"

using namespace kfftest;

/* The blocks of the sections are read from the read cache of the file (64 KB) when they are
 * entirely in it, through Kff_file::read() otherwise. These tests read files of several cache
 * sizes made of blocks of random sizes (some larger than the cache) with each read function, and
 * check the blocks and the position in the file after each of them. */

/** Sizes of blocks: mostly small, sometimes very large */
static uint64_t random_nb_kmers(Random & rnd, uint64_t max) {
	const uint64_t r = rnd.between(0, 99);
	if (r == 0)
		return max;
	if (r < 10)
		return rnd.between(1, std::min<uint64_t>(max, 3000));
	return rnd.between(1, std::min<uint64_t>(max, 30));
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

TEST_CASE("block reads: raw blocks around the limits of the read cache") {
	TempDir dir;
	Random rnd(400);
	int file_idx = 0;
	for (uint64_t k : {3, 31, 64})
		for (uint64_t data_size : {0, 1, 4}) {
			const uint64_t max = 70000;
			std::string name = dir.file("raw" + std::to_string(file_idx++) + ".kff");
			std::vector<Block> blocks;
			for (int b = 0; b < 1500; b++)
				blocks.push_back(random_block(rnd, k, random_nb_kmers(rnd, max), data_size));
			{
				Kff_file file(name, "w");
				file.set_indexation(false);
				write_variables(file, k, 0, max, data_size);
				Section_Raw sr(&file);
				for (const Block & block : blocks) {
					Bytes seq = encode(block.sequence);
					Bytes data = block.data;
					sr.write_compacted_sequence(seq.data(), block.sequence.size(), data.data());
				}
				sr.close();
				file.close();
			}

			for (int function = 0; function < 3; function++) {
				Kff_file file(name, "r");
				Section_GV sgv(&file);
				Section_Raw sr(&file);
				Bytes seq((k + max + 2) / 4 + 1);
				Bytes data(max * data_size + 1);
				Bytes seq_data((k + max + 2) / 4 + max * data_size + 1);
				unsigned long position = file.tellp();
				for (const Block & block : blocks) {
					const uint64_t seq_bytes = (block.sequence.size() + 3) / 4;
					const uint64_t counter_bytes = bytes_for_max(max);
					if (function == 2) {
						sr.jump_sequence();
					} else {
						uint64_t nb_kmers = 0;
						if (function == 0) {
							nb_kmers = sr.read_compacted_sequence(seq.data(), data.data());
							CHECK_EQ(slice(data, 0, block.data.size()), block.data);
						} else {
							nb_kmers = sr.read_compacted_sequence(seq_data.data());
							seq = slice(seq_data, 0, seq_bytes);
							CHECK_EQ(slice(seq_data, seq_bytes, seq_bytes + block.data.size()), block.data);
						}
						CHECK_EQ(nb_kmers, block.sequence.size() - k + 1);
						CHECK_EQ(decode(seq, block.sequence.size()), block.sequence);
					}
					position += counter_bytes + seq_bytes + block.data.size();
					CHECK_EQ(file.tellp(), position);
				}
				CHECK_EQ(sr.remaining_blocks, (uint64_t)0);
				sr.close();
			}
		}
}

TEST_CASE("block reads: minimizer blocks around the limits of the read cache") {
	TempDir dir;
	Random rnd(401);
	int file_idx = 0;
	for (uint64_t m : {1, 9})
		for (uint64_t k : {m, (uint64_t)31, (uint64_t)70})
			for (uint64_t data_size : {0, 3}) {
				const uint64_t max = 70000;
				std::string name = dir.file("mini" + std::to_string(file_idx++) + ".kff");
				std::string minimizer = rnd.sequence(m);
				std::vector<Block> blocks;
				for (int b = 0; b < 1200; b++)
					blocks.push_back(random_block(rnd, k, random_nb_kmers(rnd, max), data_size, minimizer));
				{
					Kff_file file(name, "w");
					file.set_indexation(false);
					write_variables(file, k, m, max, data_size);
					Section_Minimizer sm(&file);
					Bytes mini = encode(minimizer);
					sm.write_minimizer(mini.data());
					for (const Block & block : blocks) {
						Bytes seq = encode(block.sequence);
						Bytes data = block.data;
						sm.write_compacted_sequence(seq.data(), block.sequence.size(), block.mini_pos, data.data());
					}
					sm.close();
					file.close();
				}

				for (int function = 0; function < 4; function++) {
					Kff_file file(name, "r");
					Section_GV sgv(&file);
					Section_Minimizer sm(&file);
					Bytes seq((k + max + 2) / 4 + 1);
					Bytes data(max * data_size + 1);
					Bytes seq_data((k + max + 2) / 4 + max * data_size + 8);
					unsigned long position = file.tellp();
					for (const Block & block : blocks) {
						const std::string without = block.sequence.substr(0, block.mini_pos) + block.sequence.substr(block.mini_pos + m);
						const uint64_t seq_bytes = (without.size() + 3) / 4;
						uint64_t nb_kmers = block.sequence.size() - k + 1;
						uint64_t mini_pos = 0;
						if (function == 0) {
							CHECK_EQ(sm.read_compacted_sequence_without_mini(seq.data(), data.data(), mini_pos), nb_kmers);
							CHECK_EQ(decode(seq, without.size()), without);
							CHECK_EQ(slice(data, 0, block.data.size()), block.data);
							CHECK_EQ(mini_pos, block.mini_pos);
						} else if (function == 1) {
							CHECK_EQ(sm.read_compacted_sequence_without_mini(seq_data.data(), mini_pos), nb_kmers);
							CHECK_EQ(decode(seq_data, without.size()), without);
							CHECK_EQ(slice(seq_data, seq_bytes, seq_bytes + block.data.size()), block.data);
							CHECK_EQ(mini_pos, block.mini_pos);
						} else if (function == 2) {
							CHECK_EQ(sm.read_compacted_sequence(seq.data(), data.data()), nb_kmers);
							CHECK_EQ(decode(seq, block.sequence.size()), block.sequence);
							CHECK_EQ(slice(data, 0, block.data.size()), block.data);
						} else {
							sm.jump_sequence();
						}
						position += bytes_for_max(max) + mini_pos_bytes(k, max) + seq_bytes + block.data.size();
						CHECK_EQ(file.tellp(), position);
					}
					CHECK_EQ(sm.remaining_blocks, (uint64_t)0);
					sm.close();
				}
			}
}
