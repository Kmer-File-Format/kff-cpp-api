#include <map>

#include "framework.hpp"
#include "helpers.hpp"

using namespace kfftest;

/** Writes v, r, v, m, r sections and returns the beginning of each section */
static std::map<long, char> write_indexed_file(const std::string & name, bool indexed) {
	std::map<long, char> sections;
	Random rnd(20);
	Kff_file file(name, "w");
	file.set_indexation(indexed);

	Section_GV v1(&file);
	sections[v1.beginning] = 'v';
	v1.write_var("k", 15);
	v1.write_var("max", 10);
	v1.write_var("data_size", 1);
	v1.close();

	Section_Raw r1(&file);
	sections[r1.beginning] = 'r';
	for (int i = 0; i < 5; i++) {
		Block b = random_block(rnd, 15, 1 + i, 1);
		Bytes seq = encode(b.sequence);
		r1.write_compacted_sequence(seq.data(), b.sequence.size(), b.data.data());
	}
	r1.close();

	Section_GV v2(&file);
	sections[v2.beginning] = 'v';
	v2.write_var("k", 15);
	v2.write_var("m", 6);
	v2.write_var("max", 10);
	v2.write_var("data_size", 1);
	v2.close();

	Section_Minimizer m1(&file);
	sections[m1.beginning] = 'm';
	Bytes mini = encode("ACGTAC");
	m1.write_minimizer(mini.data());
	for (int i = 0; i < 5; i++) {
		Block b = random_block(rnd, 15, 1 + i, 1, "ACGTAC");
		Bytes seq = encode(b.sequence);
		m1.write_compacted_sequence(seq.data(), b.sequence.size(), b.mini_pos, b.data.data());
	}
	m1.close();

	Section_Raw r2(&file);
	sections[r2.beginning] = 'r';
	r2.close();

	file.close();
	return sections;
}

TEST_CASE("index: the footer index lists every section") {
	TempDir dir;
	std::string name = dir.file("indexed.kff");
	std::map<long, char> sections = write_indexed_file(name, true);
	Bytes bytes = read_file(name);

	Kff_file file(name, "r");
	CHECK(file.footer != nullptr);
	CHECK_EQ(file.index.size(), (size_t)1);
	const Section_Index * si = file.index[0];
	CHECK_EQ(si->next_index, (int64_t)0);
	CHECK_EQ(si->index.size(), sections.size());

	// Positions are relative to the end of the index section
	const long index_position = file.footer->vars.at("first_index");
	const long index_end = index_position + 17 + 9 * sections.size();
	CHECK_EQ(bytes[index_position], (uint8_t)'i');
	std::map<long, char> absolute;
	for (const auto & entry : si->index) {
		absolute[index_end + entry.first] = entry.second;
		CHECK(entry.first < 0);  // all the sections are before the index
	}
	CHECK(absolute == sections);
	for (const auto & entry : sections)
		CHECK_EQ(bytes[entry.first], (uint8_t)entry.second);

	// Exact bytes of the index section (entries sorted by relative position)
	Bytes expected = {'i'};
	expected += big_endian(sections.size(), 8);
	for (const auto & entry : sections) {
		expected += Bytes{(uint8_t)entry.second};
		expected += big_endian(entry.first - index_end, 8);
	}
	expected += big_endian(0, 8);
	CHECK_EQ(slice(bytes, index_position, index_end), expected);
	// Then the footer and the signature
	Bytes footer = variables_bytes({{"first_index", (uint64_t)index_position}, {"footer_size", 49}});
	footer += Bytes{'K', 'F', 'F'};
	CHECK_EQ(slice(bytes, index_end, bytes.size()), footer);
}

TEST_CASE("index: no footer when the indexation is disabled") {
	TempDir dir;
	std::string name = dir.file("not_indexed.kff");
	std::map<long, char> sections = write_indexed_file(name, false);
	Bytes bytes = read_file(name);
	// The last section (an empty raw section) is directly followed by the signature
	CHECK_EQ(bytes.size(), (size_t)(sections.rbegin()->first + 9 + 3));

	Kff_file file(name, "r");
	CHECK(file.footer == nullptr);
	CHECK(file.index.empty());
	// The reading starts at the first section
	CHECK_EQ(file.read_section_type(), 'v');
	CHECK_EQ(file.tellp(), (unsigned long)first_section_position());
}

TEST_CASE("index: chained index sections are all loaded") {
	TempDir dir;
	std::string name = dir.file("chained.kff");
	{
		Kff_file file(name, "w");
		file.set_indexation(false);
		// First index at the beginning of the file, pointing to a second one
		Section_Index first(&file);
		first.register_section('v', 0);
		// first index: 1 + 8 + 9 + 8 = 26 bytes, then a 'v' section of 9 bytes
		first.set_next_index(9);
		first.close();
		Section_GV sgv(&file);
		sgv.close();
		// Second index at 47, its end is at 47 + 26 = 73 and the 'v' section at 38
		Section_Index second(&file);
		second.register_section('v', 38 - 73);
		second.close();
		file.close();
	}
	Kff_file file(name, "r");
	CHECK_EQ(file.index.size(), (size_t)2);
	CHECK_EQ(file.index[0]->next_index, (int64_t)9);
	CHECK_EQ(file.index[0]->index.at(0), 'v');
	CHECK_EQ(file.index[1]->next_index, (int64_t)0);
	CHECK_EQ(file.index[1]->index.at(38 - 73), 'v');
	// Reading still starts at the first section
	CHECK_EQ(file.read_section_type(), 'i');
}

TEST_CASE("index: sections can be built from their type") {
	TempDir dir;
	std::string name = dir.file("builder.kff");
	std::map<long, char> sections = write_indexed_file(name, true);
	Kff_file file(name, "r");
	std::string types;
	while (file.tellp() < file.end_position) {
		Section * section = SectionBuilder::build(&file);
		if (Section_Raw * sr = dynamic_cast<Section_Raw *>(section)) {
			types += 'r';
			sr->close();
		} else if (Section_Minimizer * sm = dynamic_cast<Section_Minimizer *>(section)) {
			types += 'm';
			sm->close();
		} else if (Section_GV * sgv = dynamic_cast<Section_GV *>(section)) {
			types += 'v';
			sgv->close();
		} else if (Section_Index * si = dynamic_cast<Section_Index *>(section)) {
			types += 'i';
			si->close();
		}
		delete section;
	}
	CHECK_EQ(types, std::string("vrvmriv"));
	(void)sections;
}

TEST_CASE("index: entries registered out of order or twice (the last type is kept)") {
	TempDir dir;
	std::string name = dir.file("unordered.kff");
	{
		Kff_file file(name, "w");
		file.set_indexation(false);
		Section_Index si(&file);
		si.register_section('m', 10);
		si.register_section('r', -5);
		si.register_section('v', 10);
		si.register_section('r', 30);
		si.register_section('v', -20);
		si.register_section('m', 30);
		CHECK(si.index == (std::map<int64_t, char>{{-20, 'v'}, {-5, 'r'}, {10, 'v'}, {30, 'm'}}));
		si.close();
		file.close();
	}
	Bytes bytes = read_file(name);
	Bytes expected = {'i'};
	expected += big_endian(4, 8);
	for (auto entry : std::vector<std::pair<int64_t, char>>{{-20, 'v'}, {-5, 'r'}, {10, 'v'}, {30, 'm'}}) {
		expected += Bytes{(uint8_t)entry.second};
		expected += big_endian(entry.first, 8);
	}
	expected += big_endian(0, 8);
	CHECK_EQ(slice(bytes, first_section_position(), bytes.size() - 3), expected);

	// Positions of the sections of a file registered before the last one or twice (the file
	// starts with 8 bytes of header)
	Kff_file file(dir.file("positions.kff"), "w");
	Bytes zeros(100, 0);
	file.write(zeros.data(), zeros.size());
	file.register_position('r');
	file.register_position('v');
	file.jump_to(40);
	file.register_position('r');
	file.register_position('m');
	file.jump_to(100);
	file.register_position('r');
	CHECK(file.section_positions == (std::map<long, char>{{40, 'm'}, {100, 'r'}, {108, 'v'}}));
	file.close(false);
}

TEST_CASE("index: footer of a file of many sections") {
	TempDir dir;
	std::string name = dir.file("many.kff");
	const int nb_sections = 10000;
	std::map<long, char> sections;
	{
		Kff_file file(name, "w");
		Section_GV sgv(&file);
		sections[sgv.beginning] = 'v';
		sgv.write_var("k", 10);
		sgv.write_var("m", 5);
		sgv.write_var("max", 3);
		sgv.write_var("data_size", 0);
		sgv.close();
		Random rnd(30);
		for (int s = 0; s < nb_sections; s++) {
			if (s % 3 == 0) {
				Section_Raw sr(&file);
				sections[sr.beginning] = 'r';
				Bytes seq = encode(rnd.sequence(11));
				sr.write_compacted_sequence(seq.data(), 11, nullptr);
				sr.close();
			} else {
				Section_Minimizer sm(&file);
				sections[sm.beginning] = 'm';
				Bytes mini = encode(rnd.sequence(5));
				sm.write_minimizer(mini.data());
				Bytes seq = encode(rnd.sequence(6));
				sm.write_compacted_sequence_without_mini(seq.data(), 6, 2, nullptr);
				sm.close();
			}
		}
		file.close();
	}
	Bytes bytes = read_file(name);
	Kff_file file(name, "r");
	const long index_position = file.footer->vars.at("first_index");
	const long index_end = index_position + 17 + 9 * sections.size();
	Bytes expected = {'i'};
	expected += big_endian(sections.size(), 8);
	for (const auto & entry : sections) {
		expected += Bytes{(uint8_t)entry.second};
		expected += big_endian(entry.first - index_end, 8);
	}
	expected += big_endian(0, 8);
	CHECK_EQ(slice(bytes, index_position, index_end), expected);
	CHECK_EQ(file.index[0]->index.size(), sections.size());
}
