#include <algorithm>
#include <stdexcept>

#include "framework.hpp"
#include "helpers.hpp"

using namespace kfftest;

// Size of the automatic footer: index section without entries + 'v' footer section
static const size_t EMPTY_INDEX_SIZE = 1 + 8 + 8;
static const size_t FOOTER_SIZE = 1 + 8 + (12 + 8) + (12 + 8);

TEST_CASE("header: default header of a new file") {
	TempDir dir;
	std::string name = dir.file("empty.kff");
	{
		Kff_file file(name, "w");
		file.close();
	}
	Bytes bytes = read_file(name);
	// Signature, version, default encoding (A=0 C=1 G=3 T=2), flags, empty metadata
	Bytes expected = {'K', 'F', 'F', KFF_FORMAT_VERSION_MAJOR, KFF_FORMAT_VERSION_MINOR, 0x1e, 0, 0, 0, 0, 0, 0};
	// Index with 0 sections
	expected += Bytes{'i'};
	expected += big_endian(0, 8);
	expected += big_endian(0, 8);
	// Footer
	expected += variables_bytes({{"first_index", 12}, {"footer_size", FOOTER_SIZE}});
	expected += Bytes{'K', 'F', 'F'};
	CHECK_EQ(bytes, expected);
	CHECK_EQ(bytes.size(), 12 + EMPTY_INDEX_SIZE + FOOTER_SIZE + 3);
}

TEST_CASE("header: encodings, flags and metadata round trip") {
	TempDir dir;
	uint8_t encoding[4] = {0, 1, 2, 3};
	const Bytes metadata = {'h', 'e', 'l', 'l', 'o', 0, 0xff, 'w', 'o', 'r', 'l', 'd'};
	int permutation = 0;
	do {
		bool unique = permutation % 2 == 0;
		bool canonical = permutation % 3 == 0;
		std::string name = dir.file("header" + std::to_string(permutation) + ".kff");
		{
			Kff_file file(name, "w");
			file.write_encoding(encoding[0], encoding[1], encoding[2], encoding[3]);
			file.set_uniqueness(unique);
			file.set_canonicity(canonical);
			file.write_metadata(metadata.size(), metadata.data());
			file.close();
		}
		Bytes bytes = read_file(name);
		CHECK_EQ(slice(bytes, 0, 12 + metadata.size()), header_bytes(encoding, unique, canonical, metadata));

		Kff_file file(name, "r");
		for (int i = 0; i < 4; i++)
			CHECK_EQ(file.encoding[i], encoding[i]);
		CHECK_EQ(file.uniqueness, unique);
		CHECK_EQ(file.canonicity, canonical);
		CHECK_EQ(file.major_version, (uint8_t)KFF_FORMAT_VERSION_MAJOR);
		CHECK_EQ(file.minor_version, (uint8_t)KFF_FORMAT_VERSION_MINOR);
		CHECK_EQ(file.metadata_size, (uint32_t)metadata.size());
		Bytes read(file.metadata_size);
		file.read_metadata(read.data());
		CHECK_EQ(read, metadata);
		permutation += 1;
	} while (std::next_permutation(encoding, encoding + 4));
	CHECK_EQ(permutation, 24);
}

TEST_CASE("header: write_encoding(uint8_t*) and only the 2 lower bits are used") {
	TempDir dir;
	std::string name = dir.file("enc.kff");
	{
		Kff_file file(name, "w");
		uint8_t encoding[4] = {0xff, 0xfe, 0xfc, 0xfd};  // = 3 2 0 1
		file.write_encoding(encoding);
		file.close();
	}
	Kff_file file(name, "r");
	CHECK_EQ(file.encoding[0], (uint8_t)3);
	CHECK_EQ(file.encoding[1], (uint8_t)2);
	CHECK_EQ(file.encoding[2], (uint8_t)0);
	CHECK_EQ(file.encoding[3], (uint8_t)1);
}

TEST_CASE("header: an invalid encoding is rejected when written") {
	TempDir dir;
	Kff_file file(dir.file("bad.kff"), "w");
	CHECK_THROWS_AS(file.write_encoding(0, 0, 1, 2), std::invalid_argument);
	CHECK_THROWS_AS(file.write_encoding(0, 1, 0, 2), std::invalid_argument);
	CHECK_THROWS_AS(file.write_encoding(0, 1, 2, 0), std::invalid_argument);
	CHECK_THROWS_AS(file.write_encoding(0, 1, 1, 2), std::invalid_argument);
	CHECK_THROWS_AS(file.write_encoding(0, 1, 2, 1), std::invalid_argument);
	CHECK_THROWS_AS(file.write_encoding(0, 1, 2, 2), std::invalid_argument);
	file.write_encoding(3, 2, 1, 0);
	file.close();
}

/** Write a minimal valid file and return its bytes */
static Bytes valid_file_bytes(const TempDir & dir) {
	std::string name = dir.file("valid.kff");
	{
		Kff_file file(name, "w");
		file.close();
	}
	return read_file(name);
}

TEST_CASE("header: an invalid encoding is rejected when read") {
	TempDir dir;
	Bytes bytes = valid_file_bytes(dir);
	bytes[5] = 0b00011100;  // A=0 C=1 G=3 T=0
	write_file(dir.file("bad.kff"), bytes);
	CHECK_THROWS_AS(Kff_file(dir.file("bad.kff"), "r"), const char *);
}

TEST_CASE("header: missing signatures are detected") {
	TempDir dir;
	Bytes bytes = valid_file_bytes(dir);

	Bytes start = bytes;
	start[1] = 'X';
	write_file(dir.file("start.kff"), start);
	CHECK_THROWS_AS(Kff_file(dir.file("start.kff"), "r"), const char *);

	Bytes end = bytes;
	end[end.size() - 1] = 'X';
	write_file(dir.file("end.kff"), end);
	CHECK_THROWS_AS(Kff_file(dir.file("end.kff"), "r"), const char *);
}

TEST_CASE("header: the files are written in the version 1.0 of the format") {
	// The version of the format only changes with the format, not with the version of the library
	CHECK_EQ(KFF_FORMAT_VERSION_MAJOR, 1);
	CHECK_EQ(KFF_FORMAT_VERSION_MINOR, 0);
	// Former names of the version of the format
	CHECK_EQ(KFF_VERSION_MAJOR, KFF_FORMAT_VERSION_MAJOR);
	CHECK_EQ(KFF_VERSION_MINOR, KFF_FORMAT_VERSION_MINOR);
	TempDir dir;
	std::string name = dir.file("version.kff");
	{
		Kff_file file(name, "w");
		file.close();
	}
	Bytes bytes = read_file(name);
	CHECK_EQ(bytes[3], (uint8_t)1);
	CHECK_EQ(bytes[4], (uint8_t)0);
}

TEST_CASE("header: files from a newer format version are rejected") {
	TempDir dir;
	Bytes bytes = valid_file_bytes(dir);

	Bytes major = bytes;
	major[3] = KFF_FORMAT_VERSION_MAJOR + 1;
	write_file(dir.file("major.kff"), major);
	CHECK_THROWS_AS(Kff_file(dir.file("major.kff"), "r"), const char *);

	Bytes minor = bytes;
	minor[4] = KFF_FORMAT_VERSION_MINOR + 1;
	write_file(dir.file("minor.kff"), minor);
	CHECK_THROWS_AS(Kff_file(dir.file("minor.kff"), "r"), const char *);
}

TEST_CASE("header: a missing file cannot be opened") {
	TempDir dir;
	CHECK_THROWS_AS(Kff_file(dir.file("absent.kff"), "r"), std::runtime_error);
}

TEST_CASE("header: metadata is skipped when not read") {
	TempDir dir;
	std::string name = dir.file("meta.kff");
	const Bytes metadata(1000, 'x');
	{
		Kff_file file(name, "w");
		file.write_metadata(metadata.size(), metadata.data());
		Section_GV sgv(&file);
		sgv.write_var("k", 12);
		sgv.close();
		file.close();
	}
	Kff_file file(name, "r");
	CHECK_EQ(file.read_section_type(), 'v');
	Section_GV sgv(&file);
	CHECK_EQ(sgv.vars.at("k"), (uint64_t)12);
}

TEST_CASE("header: smallest valid file, without section nor footer") {
	TempDir dir;
	Bytes bytes = header_bytes();
	bytes += Bytes{'K', 'F', 'F'};
	CHECK_EQ(bytes.size(), (size_t)15);
	write_file(dir.file("small.kff"), bytes);

	Kff_file file(dir.file("small.kff"), "r");
	CHECK(file.footer == nullptr);
	CHECK_EQ(file.tellp(), (unsigned long)12);
	CHECK_EQ(file.end_position, (unsigned long)12);

	Kff_reader reader(dir.file("small.kff"));
	CHECK(!reader.has_next());
}
