#include "framework.hpp"
#include "helpers.hpp"

using namespace kfftest;

struct Kmer {
	std::string sequence;
	Bytes data;
	uint64_t k;
};

/** File with raw and minimizer sections that change k, max and data_size */
static std::vector<Kmer> write_mixed_file(const std::string & name, const uint8_t encoding[4]) {
	std::vector<Kmer> kmers;
	Random rnd(30);
	Kff_file file(name, "w");
	file.write_encoding(encoding[0], encoding[1], encoding[2], encoding[3]);

	auto add_kmers = [&](const Block & block, uint64_t k, uint64_t data_size) {
		for (size_t i = 0; i + k <= block.sequence.size(); i++)
			kmers.push_back({block.sequence.substr(i, k),
			                 slice(block.data, i * data_size, (i + 1) * data_size), k});
	};

	struct Params { uint64_t k, m, max, data_size; char type; };
	const Params params[] = {
		{5, 0, 4, 1, 'r'}, {7, 0, 3, 0, 'r'}, {6, 3, 5, 2, 'm'}, {33, 0, 1, 3, 'r'},
		{6, 3, 5, 2, 'r'}, {31, 11, 40, 1, 'm'}, {64, 0, 9, 1, 'r'}};
	for (const Params & p : params) {
		Section_GV sgv(&file);
		sgv.write_var("k", p.k);
		sgv.write_var("max", p.max);
		sgv.write_var("data_size", p.data_size);
		if (p.type == 'm')
			sgv.write_var("m", p.m);
		sgv.close();

		if (p.type == 'r') {
			Section_Raw sr(&file);
			for (int i = 0; i < 25; i++) {
				Block block = random_block(rnd, p.k, rnd.between(1, p.max), p.data_size);
				Bytes seq = encode(block.sequence, encoding);
				sr.write_compacted_sequence(seq.data(), block.sequence.size(), block.data.data());
				add_kmers(block, p.k, p.data_size);
			}
			sr.close();
			// An empty section that must be skipped
			Section_Raw empty(&file);
			empty.close();
		} else {
			Section_Minimizer sm(&file);
			std::string minimizer = rnd.sequence(p.m);
			Bytes mini = encode(minimizer, encoding);
			sm.write_minimizer(mini.data());
			for (int i = 0; i < 25; i++) {
				Block block = random_block(rnd, p.k, rnd.between(1, p.max), p.data_size, minimizer);
				Bytes seq = encode(block.sequence, encoding);
				sm.write_compacted_sequence(seq.data(), block.sequence.size(), block.mini_pos, block.data.data());
				add_kmers(block, p.k, p.data_size);
			}
			sm.close();
		}
	}
	file.close();
	return kmers;
}

TEST_CASE("reader: next_kmer enumerates every kmer and its data") {
	TempDir dir;
	const uint8_t encoding[4] = {2, 0, 3, 1};
	std::string name = dir.file("mixed.kff");
	std::vector<Kmer> expected = write_mixed_file(name, encoding);

	Kff_reader reader(name);
	CHECK_EQ(reader.get_encoding()[0], encoding[0]);
	CHECK_EQ(reader.get_encoding()[3], encoding[3]);
	uint8_t * kmer = nullptr;
	uint8_t * data = nullptr;
	size_t idx = 0;
	while (reader.has_next()) {
		CHECK(reader.next_kmer(kmer, data));
		CHECK(idx < expected.size());
		CHECK_EQ(reader.k, expected[idx].k);
		CHECK_EQ(decode(kmer, reader.k, encoding), expected[idx].sequence);
		CHECK_EQ(Bytes(data, data + reader.data_size), expected[idx].data);
		idx += 1;
	}
	CHECK_EQ(idx, expected.size());
	CHECK(!reader.next_kmer(kmer, data));
	CHECK(kmer == nullptr);
	CHECK(data == nullptr);
}

TEST_CASE("reader: next_block returns whole blocks") {
	TempDir dir;
	std::string name = dir.file("blocks.kff");
	Random rnd(31);
	std::vector<Block> blocks;
	{
		Kff_file file(name, "w");
		Section_GV sgv(&file);
		sgv.write_var("k", 13);
		sgv.write_var("max", 50);
		sgv.write_var("data_size", 2);
		sgv.close();
		Section_Raw sr(&file);
		for (int i = 0; i < 40; i++) {
			blocks.push_back(random_block(rnd, 13, rnd.between(1, 50), 2));
			Bytes seq = encode(blocks.back().sequence);
			sr.write_compacted_sequence(seq.data(), blocks.back().sequence.size(), blocks.back().data.data());
		}
		sr.close();
		file.close();
	}
	Kff_reader reader(name);
	CHECK_EQ(reader.get_var("max"), (uint64_t)50);
	CHECK_EQ(reader.get_var("data_size"), (uint64_t)2);
	Bytes seq_buffer(20), data_buffer(100);
	uint8_t * seq = seq_buffer.data();
	uint8_t * data = data_buffer.data();
	for (const Block & block : blocks) {
		uint64_t nb_kmers = reader.next_block(seq, data);
		CHECK_EQ(nb_kmers, block.sequence.size() - 12);
		CHECK_EQ(decode(seq, block.sequence.size()), block.sequence);
		CHECK_EQ(Bytes(data, data + block.data.size()), block.data);
	}
	CHECK(!reader.has_next());
	CHECK_EQ(reader.next_block(seq, data), (uint64_t)0);
	CHECK(seq == nullptr);
}

TEST_CASE("reader: file without kmer") {
	TempDir dir;
	std::string name = dir.file("nokmer.kff");
	{
		Kff_file file(name, "w");
		Section_GV sgv(&file);
		sgv.write_var("k", 13);
		sgv.close();
		file.close();
	}
	Kff_reader reader(name);
	CHECK(!reader.has_next());
	uint8_t * kmer = nullptr;
	uint8_t * data = nullptr;
	CHECK(!reader.next_kmer(kmer, data));
}
