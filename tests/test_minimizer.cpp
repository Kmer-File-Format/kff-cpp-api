#include "framework.hpp"
#include "helpers.hpp"

using namespace kfftest;

struct MiniFile {
	long beginning;
	std::string minimizer;
	std::vector<Block> blocks;
};

/** Write a file with one 'v' and one 'm' section. Odd blocks are written with
 *  write_compacted_sequence (minimizer removed by the API), even blocks with
 *  write_compacted_sequence_without_mini. */
static MiniFile write_mini_file(const std::string & name, Random & rnd, uint64_t k, uint64_t m,
                                uint64_t max, uint64_t data_size, size_t nb_blocks) {
	MiniFile result;
	result.minimizer = rnd.sequence(m);
	Kff_file file(name, "w");
	file.set_indexation(false);
	Section_GV sgv(&file);
	sgv.write_var("k", k);
	sgv.write_var("m", m);
	sgv.write_var("max", max);
	sgv.write_var("data_size", data_size);
	sgv.close();

	Section_Minimizer sm(&file);
	result.beginning = sm.beginning;
	Bytes mini = encode(result.minimizer);
	sm.write_minimizer(mini.data());
	for (size_t i = 0; i < nb_blocks; i++) {
		uint64_t nb_kmers = i == 0 ? max : rnd.between(1, max);
		Block block = random_block(rnd, k, nb_kmers, data_size, result.minimizer);
		Bytes data = block.data;
		if (i % 2 == 1) {
			Bytes seq = encode(block.sequence);
			sm.write_compacted_sequence(seq.data(), block.sequence.size(), block.mini_pos, data.data());
		} else {
			std::string without = block.sequence.substr(0, block.mini_pos) + block.sequence.substr(block.mini_pos + m);
			Bytes seq = encode(without);
			sm.write_compacted_sequence_without_mini(seq.data(), without.size(), block.mini_pos, data.data());
		}
		result.blocks.push_back(block);
	}
	CHECK_EQ(sm.nb_blocks, (uint64_t)nb_blocks);
	sm.close();
	file.close();
	return result;
}

TEST_CASE("minimizer: section bytes and round trip for many k, m, max and data sizes") {
	TempDir dir;
	Random rnd(10);
	const uint64_t ms[] = {1, 2, 3, 4, 5, 7, 8, 9, 13, 16, 17, 31};
	const uint64_t deltas[] = {0, 1, 2, 3, 4, 5, 9, 20};  // k - m
	const uint64_t maxs[] = {1, 3, 7, 200};
	const uint64_t data_sizes[] = {0, 1, 3};
	int idx = 0;
	for (uint64_t m : ms)
		for (uint64_t delta : deltas)
			for (uint64_t max : maxs)
				for (uint64_t data_size : data_sizes) {
					const uint64_t k = m + delta;
					std::string name = dir.file("mini" + std::to_string(idx++) + ".kff");
					MiniFile mf = write_mini_file(name, rnd, k, m, max, data_size, 12);

					// Exact bytes
					Bytes bytes = read_file(name);
					CHECK_EQ(slice(bytes, mf.beginning, bytes.size() - 3),
					         minimizer_section_bytes(mf.minimizer, mf.blocks, k, max));

					const size_t seq_capacity = (k + max) / 4 + 2;
					const size_t data_capacity = max * data_size + 1;

					// 1 - Read without the minimizer
					{
						Kff_file file(name, "r");
						Section_GV sgv(&file);
						Section_Minimizer sm(&file);
						CHECK_EQ(sm.k, k);
						CHECK_EQ(sm.m, (uint8_t)m);
						CHECK_EQ(decode(sm.minimizer, m), mf.minimizer);
						Bytes seq(seq_capacity), data(data_capacity);
						for (const Block & block : mf.blocks) {
							uint64_t mini_pos = 0;
							uint64_t nb_kmers = sm.read_compacted_sequence_without_mini(seq.data(), data.data(), mini_pos);
							std::string without = block.sequence.substr(0, block.mini_pos) + block.sequence.substr(block.mini_pos + m);
							CHECK_EQ(nb_kmers, block.sequence.size() - k + 1);
							CHECK_EQ(mini_pos, block.mini_pos);
							CHECK_EQ(slice(seq, 0, (without.size() + 3) / 4), encode(without));
							CHECK_EQ(slice(data, 0, block.data.size()), block.data);
						}
					}
					// 2 - Read with the minimizer inserted back
					{
						Kff_file file(name, "r");
						Section_GV sgv(&file);
						Section_Minimizer sm(&file);
						Bytes seq(seq_capacity), data(data_capacity);
						for (const Block & block : mf.blocks) {
							uint64_t nb_kmers = sm.read_compacted_sequence(seq.data(), data.data());
							CHECK_EQ(nb_kmers, block.sequence.size() - k + 1);
							// Same bytes as the reference encoding: padding bits are set to 0
							CHECK_EQ(slice(seq, 0, (block.sequence.size() + 3) / 4), encode(block.sequence));
							CHECK_EQ(slice(data, 0, block.data.size()), block.data);
						}
					}
					// 3 - Read sequence and data in the same array
					{
						Kff_file file(name, "r");
						Section_GV sgv(&file);
						Section_Minimizer sm(&file);
						Bytes buffer(seq_capacity + data_capacity);
						for (const Block & block : mf.blocks) {
							uint64_t nb_kmers = sm.read_compacted_sequence(buffer.data());
							const size_t seq_bytes = (block.sequence.size() + 3) / 4;
							CHECK_EQ(nb_kmers, block.sequence.size() - k + 1);
							CHECK_EQ(slice(buffer, 0, seq_bytes), encode(block.sequence));
							CHECK_EQ(slice(buffer, seq_bytes, seq_bytes + block.data.size()), block.data);
						}
					}
				}
}

TEST_CASE("minimizer: width of the minimizer position") {
	// ceil(log2(k + max - 1)) bits
	const uint64_t cases[][3] = {
		// k, max, expected bytes
		{1, 1, 0}, {2, 1, 1}, {3, 255, 2}, {10, 247, 1}, {10, 248, 2}, {31, 226, 1},
		{31, 227, 2}, {31, 65506, 2}, {31, 65507, 3}};
	TempDir dir;
	Random rnd(11);
	for (const auto & c : cases) {
		const uint64_t k = c[0], max = c[1];
		CHECK_EQ(mini_pos_bytes(k, max), (size_t)c[2]);
		std::string name = dir.file("width.kff");
		MiniFile mf = write_mini_file(name, rnd, k, 1, max, 0, 1);
		Bytes bytes = read_file(name);
		CHECK_EQ(slice(bytes, mf.beginning, bytes.size() - 3),
		         minimizer_section_bytes(mf.minimizer, mf.blocks, k, max));
	}
}

TEST_CASE("minimizer: jump over blocks") {
	TempDir dir;
	Random rnd(12);
	std::string name = dir.file("jump.kff");
	MiniFile mf = write_mini_file(name, rnd, 21, 7, 30, 2, 10);
	Kff_file file(name, "r");
	Section_GV sgv(&file);
	Section_Minimizer sm(&file);
	sm.jump_sequence();
	Bytes seq(20), data(100);
	uint64_t nb_kmers = sm.read_compacted_sequence(seq.data(), data.data());
	CHECK_EQ(nb_kmers, mf.blocks[1].sequence.size() - 20);
	CHECK_EQ(decode(seq, mf.blocks[1].sequence.size()), mf.blocks[1].sequence);
	sm.close();
	CHECK_EQ(file.tellp(), file.end_position);
}

TEST_CASE("minimizer: copy to another file") {
	TempDir dir;
	Random rnd(13);
	std::string source = dir.file("source.kff");
	write_mini_file(source, rnd, 25, 9, 40, 1, 30);
	std::string copy = dir.file("copy.kff");
	{
		Kff_file in(source, "r");
		Kff_file out(copy, "w");
		out.set_indexation(false);
		Section_GV sgv(&in);
		sgv.copy(&out);
		Section_Minimizer sm(&in);
		sm.copy(&out);
		out.close();
	}
	CHECK_EQ(read_file(copy), read_file(source));
}

TEST_CASE("minimizer: missing m is detected") {
	TempDir dir;
	Kff_file file(dir.file("missing.kff"), "w");
	Section_GV sgv(&file);
	sgv.write_var("k", 5);
	sgv.write_var("max", 5);
	sgv.write_var("data_size", 0);
	sgv.close();
	CHECK_THROWS_AS(Section_Minimizer section(&file), const char *);
	file.close();
}

TEST_CASE("minimizer: the README example") {
	// k=10, m=8, max=255, data_size=1, encoding A=0 C=2 G=3 T=1
	const uint8_t encoding[4] = {0, 2, 3, 1};
	TempDir dir;
	std::string name = dir.file("readme_mini.kff");
	long beginning = 0;
	{
		Kff_file file(name, "w");
		file.set_indexation(false);
		file.write_encoding(0, 2, 3, 1);
		Section_GV sgv(&file);
		sgv.write_var("k", 10);
		sgv.write_var("m", 8);
		sgv.write_var("max", 255);
		sgv.write_var("data_size", 1);
		sgv.close();
		Section_Minimizer sm(&file);
		beginning = sm.beginning;
		Bytes mini = encode("AAACTGAT", encoding);
		sm.write_minimizer(mini.data());
		Bytes s1 = encode("ACTAAACTGATT", encoding), d1 = {32, 47, 1};
		Bytes s2 = encode("AAACTGATCG", encoding), d2 = {12};
		Bytes s3 = encode("CTAAACTGATT", encoding), d3 = {1, 47};
		sm.write_compacted_sequence(s1.data(), 12, 3, d1.data());
		sm.write_compacted_sequence(s2.data(), 10, 0, d2.data());
		sm.write_compacted_sequence(s3.data(), 11, 2, d3.data());
		sm.close();
		file.close();
	}
	Bytes bytes = read_file(name);
	// Same blocks as the README, except that the minimizer position uses
	// ceil(log2(k + max - 1)) = 9 bits, i.e. 2 bytes (the README draws 1 byte).
	const Bytes expected = {'m', 0x02, 0x71, 0, 0, 0, 0, 0, 0, 0, 3,
	                        3, 0, 3, 0x25, 0x20, 0x2f, 0x01,
	                        1, 0, 0, 0x0b, 0x0c,
	                        2, 0, 2, 0x25, 0x01, 0x2f};
	CHECK_EQ(slice(bytes, beginning, bytes.size() - 3), expected);

	// The full sequences are rebuilt when reading
	Kff_file file(name, "r");
	Section_GV sgv(&file);
	Section_Minimizer sm(&file);
	Bytes seq(4), data(3);
	CHECK_EQ(sm.read_compacted_sequence(seq.data(), data.data()), (uint64_t)3);
	CHECK_EQ(decode(seq, 12, encoding), std::string("ACTAAACTGATT"));
	CHECK_EQ(sm.read_compacted_sequence(seq.data(), data.data()), (uint64_t)1);
	CHECK_EQ(decode(seq, 10, encoding), std::string("AAACTGATCG"));
	CHECK_EQ(sm.read_compacted_sequence(seq.data(), data.data()), (uint64_t)2);
	CHECK_EQ(decode(seq, 11, encoding), std::string("CTAAACTGATT"));
}
