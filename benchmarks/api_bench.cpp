#include <chrono>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include "kff_io.hpp"

/* Benchmarks of the different layers of the API on an existing kff file.
 *
 *   api_bench kmers  in.kff           high level API, kmer per kmer
 *   api_bench blocks in.kff           low level API, block per block (minimizers re-inserted)
 *   api_bench copy   in.kff out.kff   section per section copy (Section::copy)
 *   api_bench remini in.kff out.kff   rewrite minimizer sections from full sequences
 *                                     (write_compacted_sequence removes the minimizer)
 *
 * The time and a checksum of what was read are printed. The output files of
 * copy and remini can be compared between two versions of the API.
 */

using namespace std;

static double now() {
	return chrono::duration<double>(chrono::steady_clock::now().time_since_epoch()).count();
}

static uint64_t bench_kmers(const string & filename) {
	Kff_reader reader(filename);
	uint8_t * kmer;
	uint8_t * data;
	uint64_t checksum = 0;
	uint64_t nb_kmers = 0;
	while (reader.next_kmer(kmer, data)) {
		const uint64_t k = reader.k;
		// Only the nucleotides: the padding bits of the first byte are not specified
		const uint8_t mask = (k % 4 == 0) ? 0xff : (uint8_t)((1u << (2 * (k % 4))) - 1);
		checksum = checksum * 31 + (kmer[0] & mask);
		for (uint64_t b = 1; b < (k + 3) / 4; b++)
			checksum = checksum * 31 + kmer[b];
		for (uint64_t b = 0; b < reader.data_size; b++)
			checksum = checksum * 31 + data[b];
		nb_kmers += 1;
	}
	cerr << nb_kmers << " kmers" << endl;
	return checksum;
}

static uint64_t bench_blocks(const string & filename) {
	Kff_file file(filename, "r");
	vector<uint8_t> seq, data;
	uint64_t checksum = 0;
	uint64_t nb_blocks = 0;
	while (file.tellp() < file.end_position) {
		char type = file.read_section_type();
		if (type == 'v') {
			Section_GV sgv(&file);
			sgv.close();
		} else if (type == 'i') {
			Section_Index si(&file);
			si.close();
		} else if (type == 'r' or type == 'm') {
			Block_section_reader * section = Block_section_reader::construct_section(&file);
			seq.resize((section->k + section->max) / 4 + 2);
			data.resize(section->max * section->data_size + 1);
			const uint64_t blocks = section->nb_blocks;
			for (uint64_t b = 0; b < blocks; b++) {
				uint64_t n = section->read_compacted_sequence(seq.data(), data.data());
				for (uint64_t i = 0; i < (n + section->k + 2) / 4; i++)
					checksum = checksum * 31 + seq[i];
				for (uint64_t i = 0; i < n * section->data_size; i++)
					checksum = checksum * 31 + data[i];
				nb_blocks += 1;
			}
			delete section;
		} else {
			cerr << "Unknown section " << type << endl;
			exit(1);
		}
	}
	cerr << nb_blocks << " blocks" << endl;
	return checksum;
}

static void header_copy(Kff_file & in, Kff_file & out) {
	out.write_encoding(in.encoding);
	out.set_uniqueness(in.uniqueness);
	out.set_canonicity(in.canonicity);
	vector<uint8_t> metadata(in.metadata_size + 1);
	in.read_metadata(metadata.data());
	out.write_metadata(in.metadata_size, metadata.data());
}

static uint64_t bench_copy(const string & in_name, const string & out_name, bool remini) {
	Kff_file in(in_name, "r");
	Kff_file out(out_name, "w");
	header_copy(in, out);
	vector<uint8_t> seq, data;
	while (in.tellp() < in.end_position) {
		char type = in.read_section_type();
		if (type == 'v') {
			Section_GV sgv(&in);
			sgv.vars.erase("first_index");
			sgv.vars.erase("footer_size");
			sgv.copy(&out);
		} else if (type == 'i') {
			Section_Index si(&in);
			si.close();
		} else if (type == 'r') {
			Section_Raw sr(&in);
			sr.copy(&out);
		} else if (type == 'm' and not remini) {
			Section_Minimizer sm(&in);
			sm.copy(&out);
		} else if (type == 'm') {
			Section_Minimizer sm(&in);
			Section_Minimizer osm(&out);
			osm.write_minimizer(sm.minimizer);
			seq.resize((sm.k + sm.max) / 4 + 2);
			data.resize(sm.max * sm.data_size + 1);
			vector<uint8_t> without(seq.size());
			for (uint64_t b = 0; b < sm.nb_blocks; b++) {
				// Read the position of the minimizer, then the full sequence
				uint64_t mini_pos;
				uint64_t position = in.tellp();
				sm.read_compacted_sequence_without_mini(without.data(), data.data(), mini_pos);
				in.jump_to(position);
				sm.remaining_blocks += 1;
				uint64_t n = sm.read_compacted_sequence(seq.data(), data.data());
				osm.write_compacted_sequence(seq.data(), n + sm.k - 1, mini_pos, data.data());
			}
			osm.close();
		}
	}
	out.close();
	return 0;
}

int main(int argc, char * argv[]) {
	if (argc < 3) {
		cerr << "usage: api_bench kmers|blocks in.kff  or  api_bench copy|remini in.kff out.kff" << endl;
		return 1;
	}
	string mode = argv[1];
	double start = now();
	uint64_t checksum = 0;
	if (mode == "kmers")
		checksum = bench_kmers(argv[2]);
	else if (mode == "blocks")
		checksum = bench_blocks(argv[2]);
	else if ((mode == "copy" or mode == "remini") and argc == 4)
		checksum = bench_copy(argv[2], argv[3], mode == "remini");
	else {
		cerr << "unknown mode " << mode << endl;
		return 1;
	}
	cout << mode << "\t" << (now() - start) << " s\tchecksum " << checksum << endl;
	return 0;
}
