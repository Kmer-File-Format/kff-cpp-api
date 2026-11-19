#include "framework.hpp"
#include "helpers.hpp"

using namespace kfftest;

/** Write a file made of one 'v' section and one 'r' section (no index) */
static long write_raw_file(const std::string & name, uint64_t k, uint64_t max, uint64_t data_size,
                           const std::vector<Block> & blocks) {
	Kff_file file(name, "w");
	file.set_indexation(false);
	Section_GV sgv(&file);
	sgv.write_var("k", k);
	sgv.write_var("max", max);
	sgv.write_var("data_size", data_size);
	sgv.close();

	Section_Raw sr(&file);
	long beginning = sr.beginning;
	for (const Block & block : blocks) {
		Bytes seq = encode(block.sequence);
		Bytes data = block.data;
		sr.write_compacted_sequence(seq.data(), block.sequence.size(), data.data());
	}
	CHECK_EQ(sr.nb_blocks, (uint64_t)blocks.size());
	sr.close();
	file.close();
	return beginning;
}

static std::vector<Block> random_blocks(Random & rnd, uint64_t k, uint64_t max, uint64_t data_size, size_t nb) {
	std::vector<Block> blocks;
	for (size_t i = 0; i < nb; i++)
		blocks.push_back(random_block(rnd, k, rnd.between(1, max), data_size));
	return blocks;
}

TEST_CASE("raw: section bytes and round trip for many k, max and data sizes") {
	TempDir dir;
	Random rnd(1);
	const uint64_t ks[] = {1, 2, 3, 4, 5, 7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64, 65, 127};
	// Powers of 256 are avoided: n == max does not fit in the counter field (see the README)
	const uint64_t maxs[] = {1, 2, 3, 5, 100, 255, 300};
	const uint64_t data_sizes[] = {0, 1, 2, 3, 8};
	int file_idx = 0;
	for (uint64_t k : ks)
		for (uint64_t max : maxs)
			for (uint64_t data_size : data_sizes) {
				std::string name = dir.file("raw" + std::to_string(file_idx++) + ".kff");
				std::vector<Block> blocks = random_blocks(rnd, k, max, data_size, 20);
				blocks.push_back(random_block(rnd, k, max, data_size));  // a full block
				long beginning = write_raw_file(name, k, max, data_size, blocks);

				// Exact bytes
				Bytes bytes = read_file(name);
				Bytes expected = raw_section_bytes(blocks, k, max);
				CHECK_EQ(slice(bytes, beginning, bytes.size() - 3), expected);

				// Read back
				Kff_file file(name, "r");
				Section_GV sgv(&file);
				Section_Raw sr(&file);
				CHECK_EQ(sr.nb_blocks, (uint64_t)blocks.size());
				CHECK_EQ(sr.k, k);
				CHECK_EQ(sr.max, max);
				CHECK_EQ(sr.data_size, data_size);
				Bytes seq((k + max + 2) / 4 + 1);
				Bytes data(max * data_size + 1);
				for (const Block & block : blocks) {
					uint64_t nb_kmers = sr.read_compacted_sequence(seq.data(), data.data());
					CHECK_EQ(nb_kmers, block.sequence.size() - k + 1);
					CHECK_EQ(decode(seq, block.sequence.size()), block.sequence);
					CHECK_EQ(slice(data, 0, block.data.size()), block.data);
				}
				CHECK_EQ(sr.remaining_blocks, (uint64_t)0);
				sr.close();
			}
}

TEST_CASE("raw: width of the kmer counter depends on max") {
	// ceil(log2(max)) bits, 0 bytes when max = 1
	const std::pair<uint64_t, size_t> widths[] = {
		{1, 0}, {2, 1}, {3, 1}, {255, 1}, {256, 1}, {257, 2}, {65535, 2}, {65536, 2},
		{65537, 3}, {uint64_t(1) << 24, 3}, {(uint64_t(1) << 24) + 1, 4}, {uint64_t(1) << 32, 4},
		{(uint64_t(1) << 32) + 1, 5}};
	TempDir dir;
	Random rnd(2);
	for (const auto & width : widths) {
		std::string name = dir.file("width.kff");
		std::vector<Block> blocks = {random_block(rnd, 5, 1, 1), random_block(rnd, 5, 3, 1)};
		if (width.first < 3)
			blocks.pop_back();
		long beginning = write_raw_file(name, 5, width.first, 1, blocks);
		Bytes bytes = read_file(name);
		size_t expected_size = 9;
		for (const Block & block : blocks)
			expected_size += width.second + 2 + block.data.size();
		CHECK_EQ(bytes.size() - 3 - beginning, expected_size);
		CHECK_EQ(bytes_for_max(width.first), width.second);
	}
}

TEST_CASE("raw: read_compacted_sequence(seq_data) puts the data after the sequence") {
	TempDir dir;
	Random rnd(3);
	for (uint64_t k : {3, 4, 5, 6, 21}) {
		std::string name = dir.file("seqdata.kff");
		std::vector<Block> blocks = random_blocks(rnd, k, 50, 2, 30);
		write_raw_file(name, k, 50, 2, blocks);
		Kff_file file(name, "r");
		Section_GV sgv(&file);
		Section_Raw sr(&file);
		Bytes buffer((k + 50) / 4 + 1 + 100);
		for (const Block & block : blocks) {
			uint64_t nb_kmers = sr.read_compacted_sequence(buffer.data());
			size_t seq_bytes = (block.sequence.size() + 3) / 4;
			CHECK_EQ(nb_kmers, block.sequence.size() - k + 1);
			CHECK_EQ(decode(buffer, block.sequence.size()), block.sequence);
			CHECK_EQ(slice(buffer, seq_bytes, seq_bytes + block.data.size()), block.data);
		}
	}
}

TEST_CASE("raw: jump over blocks and sections") {
	TempDir dir;
	Random rnd(4);
	std::string name = dir.file("jump.kff");
	std::vector<Block> blocks = random_blocks(rnd, 11, 20, 1, 10);
	write_raw_file(name, 11, 20, 1, blocks);

	// jump_sequence skips exactly one block
	{
		Kff_file file(name, "r");
		Section_GV sgv(&file);
		Section_Raw sr(&file);
		sr.jump_sequence();
		sr.jump_sequence();
		Bytes seq(10), data(20);
		uint64_t nb_kmers = sr.read_compacted_sequence(seq.data(), data.data());
		CHECK_EQ(nb_kmers, blocks[2].sequence.size() - 10);
		CHECK_EQ(decode(seq, blocks[2].sequence.size()), blocks[2].sequence);
		CHECK_EQ(sr.remaining_blocks, (uint64_t)7);
		// close() jumps over the remaining blocks
		sr.close();
		CHECK_EQ(file.tellp(), file.end_position);
	}
	// jump_section / jump_next_section
	{
		Kff_file file(name, "r");
		Section_GV sgv(&file);
		Section_Raw sr(&file);
		sr.jump_section();
		CHECK_EQ(file.tellp(), file.end_position);
	}
	{
		Kff_file file(name, "r");
		CHECK(!file.jump_next_section());  // a 'v' section is not skipped
		Section_GV sgv(&file);
		CHECK(file.jump_next_section());
		CHECK_EQ(file.tellp(), file.end_position);
		CHECK(!file.jump_next_section());
	}
}

TEST_CASE("raw: section without block") {
	TempDir dir;
	std::string name = dir.file("noblock.kff");
	long beginning = write_raw_file(name, 7, 10, 0, {});
	Bytes bytes = read_file(name);
	CHECK_EQ(slice(bytes, beginning, bytes.size() - 3), raw_section_bytes({}, 7, 10));
	Kff_file file(name, "r");
	Section_GV sgv(&file);
	Section_Raw sr(&file);
	CHECK_EQ(sr.nb_blocks, (uint64_t)0);
	sr.close();
	CHECK_EQ(file.tellp(), file.end_position);
}

TEST_CASE("raw: missing variables are detected") {
	TempDir dir;
	Kff_file file(dir.file("missing.kff"), "w");
	Section_GV sgv(&file);
	sgv.write_var("k", 5);
	sgv.write_var("max", 5);
	sgv.close();
	CHECK_THROWS_AS(Section_Raw section(&file), const char *);
	Section_GV sgv2(&file);
	sgv2.write_var("max", 5);
	sgv2.write_var("data_size", 5);
	sgv2.close();
	CHECK_THROWS_AS(Section_Raw section(&file), const char *);
	Section_GV sgv3(&file);
	sgv3.write_var("k", 5);
	sgv3.write_var("data_size", 5);
	sgv3.close();
	CHECK_THROWS_AS(Section_Raw section(&file), const char *);
	file.close();
}

TEST_CASE("raw: copy to another file") {
	TempDir dir;
	Random rnd(5);
	std::string source = dir.file("source.kff");
	std::vector<Block> blocks = random_blocks(rnd, 9, 40, 3, 50);
	long beginning = write_raw_file(source, 9, 40, 3, blocks);

	std::string copy = dir.file("copy.kff");
	{
		Kff_file in(source, "r");
		Kff_file out(copy, "w");
		out.set_indexation(false);
		Section_GV sgv(&in);
		sgv.copy(&out);
		Section_Raw sr(&in);
		sr.copy(&out);
		out.close();
	}
	CHECK_EQ(read_file(copy), read_file(source));
	(void)beginning;
}

TEST_CASE("raw: the README example") {
	// k=10, max=255, data_size=1, encoding A=0 C=2 G=3 T=1
	const uint8_t encoding[4] = {0, 2, 3, 1};
	TempDir dir;
	std::string name = dir.file("readme_raw.kff");
	long beginning = 0;
	{
		Kff_file file(name, "w");
		file.set_indexation(false);
		file.write_encoding(0, 2, 3, 1);
		Section_GV sgv(&file);
		sgv.write_var("k", 10);
		sgv.write_var("max", 255);
		sgv.write_var("data_size", 1);
		sgv.close();
		Section_Raw sr(&file);
		beginning = sr.beginning;
		Bytes s1 = encode("ACTAAACTGATT", encoding), d1 = {32, 47, 1};
		Bytes s2 = encode("AAACTGATCG", encoding), d2 = {12};
		Bytes s3 = encode("CTAAACTGATT", encoding), d3 = {1, 47};
		sr.write_compacted_sequence(s1.data(), 12, d1.data());
		sr.write_compacted_sequence(s2.data(), 10, d2.data());
		sr.write_compacted_sequence(s3.data(), 11, d3.data());
		sr.close();
		file.close();
	}
	Bytes bytes = read_file(name);
	CHECK_EQ(bytes[5], (uint8_t)0x2d);
	// Bytes given by the README
	const Bytes expected = {'r', 0, 0, 0, 0, 0, 0, 0, 3,
	                        3, 0x24, 0x09, 0xc5, 0x20, 0x2f, 0x01,
	                        1, 0x00, 0x27, 0x1b, 0x0c,
	                        2, 0x24, 0x09, 0xc5, 0x01, 0x2f};
	CHECK_EQ(slice(bytes, beginning, bytes.size() - 3), expected);
}
