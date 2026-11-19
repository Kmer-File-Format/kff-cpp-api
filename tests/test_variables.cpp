#include <limits>

#include "framework.hpp"
#include "helpers.hpp"

using namespace kfftest;

TEST_CASE("variables: section bytes, names are written in sorted order") {
	TempDir dir;
	std::string name = dir.file("vars.kff");
	long beginning = 0;
	{
		Kff_file file(name, "w");
		file.set_indexation(false);
		Section_GV sgv(&file);
		beginning = sgv.beginning;
		sgv.write_var("max", 255);
		sgv.write_var("k", 10);
		sgv.write_var("data_size", 1);
		CHECK_EQ(sgv.nb_vars, (uint64_t)3);
		sgv.close();
		file.close();
	}
	Bytes bytes = read_file(name);
	CHECK_EQ(beginning, (long)first_section_position());
	Bytes expected = header_bytes();
	expected += variables_bytes({{"data_size", 1}, {"k", 10}, {"max", 255}});
	expected += Bytes{'K', 'F', 'F'};
	CHECK_EQ(bytes, expected);
}

TEST_CASE("variables: 64 bits values and arbitrary names round trip") {
	TempDir dir;
	std::string name = dir.file("values.kff");
	const std::vector<std::pair<std::string, uint64_t>> values = {
		{"zero", 0}, {"one", 1}, {"2^32", uint64_t(1) << 32}, {"2^63", uint64_t(1) << 63},
		{"max", std::numeric_limits<uint64_t>::max()}, {"name with spaces", 42}, {"x", 0x0123456789abcdefULL}};
	{
		Kff_file file(name, "w");
		Section_GV sgv(&file);
		for (const auto & v : values)
			sgv.write_var(v.first, v.second);
		sgv.close();
		file.close();
	}
	Kff_file file(name, "r");
	CHECK_EQ(file.read_section_type(), 'v');
	Section_GV sgv(&file);
	CHECK_EQ(sgv.nb_vars, (uint64_t)values.size());
	CHECK_EQ(sgv.vars.size(), values.size());
	for (const auto & v : values) {
		CHECK_EQ(sgv.vars.at(v.first), v.second);
		CHECK_EQ(file.global_vars.at(v.first), v.second);
	}
}

TEST_CASE("variables: a new section resets the scope of the variables") {
	TempDir dir;
	std::string name = dir.file("scopes.kff");
	{
		Kff_file file(name, "w");
		Section_GV first(&file);
		first.write_var("k", 10);
		first.write_var("m", 5);
		first.close();
		CHECK_EQ(file.global_vars.size(), (size_t)2);
		Section_GV second(&file);
		CHECK_EQ(file.global_vars.size(), (size_t)0);
		second.write_var("k", 12);
		second.close();
		CHECK_EQ(file.global_vars.at("k"), (uint64_t)12);
		file.close();
	}
	Kff_file file(name, "r");
	Section_GV first(&file);
	CHECK_EQ(file.global_vars.size(), (size_t)2);
	Section_GV second(&file);
	CHECK_EQ(file.global_vars.size(), (size_t)1);
	CHECK_EQ(file.global_vars.at("k"), (uint64_t)12);
	CHECK(file.global_vars.find("m") == file.global_vars.end());
}

TEST_CASE("variables: empty section") {
	TempDir dir;
	std::string name = dir.file("empty_vars.kff");
	{
		Kff_file file(name, "w");
		file.set_indexation(false);
		Section_GV sgv(&file);
		sgv.close();
		file.close();
	}
	Bytes bytes = read_file(name);
	Bytes expected = header_bytes();
	expected += variables_bytes({});
	expected += Bytes{'K', 'F', 'F'};
	CHECK_EQ(bytes, expected);

	Kff_file file(name, "r");
	Section_GV sgv(&file);
	CHECK_EQ(sgv.nb_vars, (uint64_t)0);
	CHECK(sgv.vars.empty());
}

TEST_CASE("variables: sections written in any order are read") {
	// The README example declares k, max then data_size (not sorted)
	TempDir dir;
	Bytes bytes = header_bytes();
	bytes += Bytes{'v'};
	bytes += big_endian(3, 8);
	bytes += c_string("k");
	bytes += big_endian(10, 8);
	bytes += c_string("max");
	bytes += big_endian(255, 8);
	bytes += c_string("data_size");
	bytes += big_endian(1, 8);
	bytes += Bytes{'K', 'F', 'F'};
	write_file(dir.file("spec.kff"), bytes);

	Kff_file file(dir.file("spec.kff"), "r");
	Section_GV sgv(&file);
	CHECK_EQ(sgv.vars.at("k"), (uint64_t)10);
	CHECK_EQ(sgv.vars.at("max"), (uint64_t)255);
	CHECK_EQ(sgv.vars.at("data_size"), (uint64_t)1);
	CHECK_EQ(file.tellp(), (unsigned long)(bytes.size() - 3));
}

TEST_CASE("variables: copy to another file") {
	TempDir dir;
	std::string source = dir.file("source.kff");
	std::string copy = dir.file("copy.kff");
	{
		Kff_file file(source, "w");
		file.set_indexation(false);
		Section_GV empty(&file);
		empty.close();
		Section_GV sgv(&file);
		sgv.write_var("k", 31);
		sgv.write_var("data_size", 2);
		sgv.close();
		file.close();
	}
	{
		Kff_file in(source, "r");
		Kff_file out(copy, "w");
		out.set_indexation(false);
		Section_GV empty(&in);
		empty.copy(&out);  // empty sections are not copied
		Section_GV sgv(&in);
		sgv.copy(&out);
		out.close();
	}
	Bytes expected = header_bytes();
	expected += variables_bytes({{"data_size", 2}, {"k", 31}});
	expected += Bytes{'K', 'F', 'F'};
	CHECK_EQ(read_file(copy), expected);
}

TEST_CASE("variables: the footer is discovered when the file is opened") {
	TempDir dir;
	std::string name = dir.file("footer.kff");
	{
		Kff_file file(name, "w");
		Section_GV sgv(&file);
		sgv.write_var("k", 21);
		sgv.close();
		file.close();
	}
	Kff_file file(name, "r");
	CHECK(file.footer != nullptr);
	CHECK_EQ(file.footer->vars.at("footer_size"), (uint64_t)49);
	// The index follows the only section of the file: 12 + 'v' section (1 + 8 + 2 + 8)
	CHECK_EQ(file.footer->vars.at("first_index"), (uint64_t)(12 + 19));
	// Opening does not move the reading position (first byte after the fixed header)
	CHECK_EQ(file.read_section_type(), 'v');
}
