#include <stdexcept>

#include "framework.hpp"
#include "helpers.hpp"

using namespace kfftest;

/** Writes a raw section large enough to go through the disk buffer several times */
static std::vector<Block> write_large_file(const std::string & name, bool tmp_close, size_t nb_blocks) {
	Random rnd(40);
	std::vector<Block> blocks;
	Kff_file file(name, "w");
	Section_GV sgv(&file);
	sgv.write_var("k", 31);
	sgv.write_var("max", 100);
	sgv.write_var("data_size", 1);
	sgv.close();
	Section_Raw sr(&file);
	for (size_t i = 0; i < nb_blocks; i++) {
		blocks.push_back(random_block(rnd, 31, rnd.between(1, 100), 1));
		Bytes seq = encode(blocks.back().sequence);
		sr.write_compacted_sequence(seq.data(), blocks.back().sequence.size(), blocks.back().data.data());
		if (tmp_close and i % 1000 == 0)
			file.tmp_close();
	}
	sr.close();
	file.close();
	return blocks;
}

TEST_CASE("file: large files are written through the buffer and read back") {
	TempDir dir;
	std::string name = dir.file("large.kff");
	std::vector<Block> blocks = write_large_file(name, false, 50000);
	Bytes bytes = read_file(name);
	CHECK(bytes.size() > 3 * (1 << 20));

	Bytes expected = header_bytes();
	expected += variables_bytes({{"data_size", 1}, {"k", 31}, {"max", 100}});
	const size_t raw_position = expected.size();
	expected += raw_section_bytes(blocks, 31, 100);
	CHECK_EQ(slice(bytes, 0, expected.size()), expected);

	Kff_file file(name, "r");
	CHECK_EQ(file.footer->vars.at("first_index"), (uint64_t)expected.size());
	Section_GV sgv(&file);
	Section_Raw sr(&file);
	CHECK_EQ(sr.beginning, (long)raw_position);
	Bytes seq(40), data(100);
	for (const Block & block : blocks) {
		sr.read_compacted_sequence(seq.data(), data.data());
		CHECK_EQ(decode(seq, block.sequence.size()), block.sequence);
		CHECK_EQ(slice(data, 0, block.data.size()), block.data);
	}
}

TEST_CASE("file: temporary closing does not change the file") {
	TempDir dir;
	write_large_file(dir.file("normal.kff"), false, 30000);
	write_large_file(dir.file("tmp_closed.kff"), true, 30000);
	CHECK(read_file(dir.file("normal.kff")) == read_file(dir.file("tmp_closed.kff")));
}

TEST_CASE("file: write_at on the disk part, the buffer part and across both") {
	TempDir dir;
	std::string name = dir.file("write_at.kff");
	Bytes content = Random(41).bytes(3 * (1 << 20) + 12345);
	{
		Kff_file file(name, "w");
		file.set_indexation(false);
		// Bytes written after the 8 bytes of the fixed header
		for (size_t pos = 0; pos < content.size(); pos += 1000)
			file.write(content.data() + pos, std::min<size_t>(1000, content.size() - pos));
		CHECK(file.file_size > 0);  // part of the file is on disk
		CHECK(file.next_free > 0);  // and part in the buffer
		const unsigned long boundary = file.file_size;
		const Bytes patch = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
		file.write_at(patch.data(), 16, 100);              // disk
		file.write_at(patch.data(), 16, boundary - 8);     // across disk and buffer
		file.write_at(patch.data(), 16, boundary + 1000);  // buffer
		for (size_t i = 0; i < 16; i++) {
			content[100 - 8 + i] = patch[i];
			content[boundary - 8 - 8 + i] = patch[i];
			content[boundary + 1000 - 8 + i] = patch[i];
		}
		CHECK_EQ(file.tellp(), (unsigned long)(8 + content.size()));
		file.close();
	}
	Bytes bytes = read_file(name);
	CHECK_EQ(bytes.size(), 8 + content.size() + 3);
	CHECK(slice(bytes, 8, 8 + content.size()) == content);
}

TEST_CASE("file: a file closed without saving can be read from memory") {
	for (size_t nb_blocks : {10, 100000}) {
		TempDir dir;
		std::string name = dir.file("memory.kff");
		Random rnd(42);
		std::vector<Block> blocks;
		{
			Kff_file file(name, "w");
			Section_GV sgv(&file);
			sgv.write_var("k", 21);
			sgv.write_var("max", 120);
			sgv.write_var("data_size", 0);
			sgv.close();
			Section_Raw sr(&file);
			for (size_t i = 0; i < nb_blocks; i++) {
				blocks.push_back(random_block(rnd, 21, rnd.between(1, 120), 0));
				Bytes seq = encode(blocks.back().sequence);
				sr.write_compacted_sequence(seq.data(), blocks.back().sequence.size(), nullptr);
			}
			sr.close();
			const bool flushed = file.file_size > 0;
			CHECK_EQ(flushed, nb_blocks > 10);

			file.close(false);
			file.open("r");
			Section_GV in_sgv(&file);
			CHECK_EQ(in_sgv.vars.at("k"), (uint64_t)21);
			Section_Raw in_sr(&file);
			CHECK_EQ(in_sr.nb_blocks, (uint64_t)nb_blocks);
			Bytes seq(40);
			for (const Block & block : blocks) {
				in_sr.read_compacted_sequence(seq.data(), nullptr);
				CHECK_EQ(decode(seq, block.sequence.size()), block.sequence);
			}
			in_sr.close();
		}
		// The partial file written on disk is removed at destruction
		CHECK(!fs::exists(name));
	}
}

TEST_CASE("file: positions, jumps and reading past the end") {
	TempDir dir;
	std::string name = dir.file("positions.kff");
	{
		Kff_file file(name, "w");
		Section_GV sgv(&file);
		sgv.write_var("k", 21);
		sgv.close();
		file.close();
	}
	Bytes bytes = read_file(name);
	Kff_file file(name, "r");
	CHECK_EQ(file.end_position, (unsigned long)(bytes.size() - 3));
	file.jump_to(3, true);
	CHECK_EQ(file.tellp(), (unsigned long)(bytes.size() - 3));
	uint8_t signature[3];
	file.read(signature, 3);
	CHECK(signature[0] == 'K' && signature[1] == 'F' && signature[2] == 'F');
	file.jump(-10);
	CHECK_EQ(file.tellp(), (unsigned long)(bytes.size() - 10));
	uint8_t buffer[16];
	CHECK_THROWS_AS(file.read(buffer, 16), std::out_of_range);
}
