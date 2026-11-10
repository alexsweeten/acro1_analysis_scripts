#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

constexpr int AA_COUNT = 20;
constexpr int HISTOGRAM_BINS = 200;
const std::string AA_ORDER = "ACDEFGHIKLMNPQRSTVWY";
const std::string BASE_ORDER = "ACGT";

struct TranslationCounts {
    std::array<uint32_t, AA_COUNT> amino_acids{};
    uint32_t stops = 0;
};

struct TestAccumulator {
    std::string frame_scheme;
    std::string ctd_reference;
    std::string normalization;
    std::string analysis_role;
    double observed_jsd = 0.0;
    uint32_t observed_sense_residues = 0;
    uint32_t observed_stop_codons = 0;
    std::vector<double> values;
};

std::vector<std::string> split_tab(const std::string& line) {
    std::vector<std::string> fields;
    size_t start = 0;
    while (true) {
        const size_t pos = line.find('\t', start);
        fields.push_back(line.substr(start, pos == std::string::npos ? pos : pos - start));
        if (pos == std::string::npos) break;
        start = pos + 1;
    }
    return fields;
}

std::unordered_map<std::string, size_t> header_index(const std::string& line) {
    const auto fields = split_tab(line);
    std::unordered_map<std::string, size_t> result;
    for (size_t index = 0; index < fields.size(); ++index) result[fields[index]] = index;
    return result;
}

int base_index(char base) {
    switch (base) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: throw std::runtime_error("Unexpected nucleotide");
    }
}

char amino_acid_for_codon(const std::string& codon) {
    static const std::unordered_map<std::string, char> table = {
        {"TTT", 'F'}, {"TTC", 'F'}, {"TTA", 'L'}, {"TTG", 'L'},
        {"TCT", 'S'}, {"TCC", 'S'}, {"TCA", 'S'}, {"TCG", 'S'},
        {"TAT", 'Y'}, {"TAC", 'Y'}, {"TAA", '*'}, {"TAG", '*'},
        {"TGT", 'C'}, {"TGC", 'C'}, {"TGA", '*'}, {"TGG", 'W'},
        {"CTT", 'L'}, {"CTC", 'L'}, {"CTA", 'L'}, {"CTG", 'L'},
        {"CCT", 'P'}, {"CCC", 'P'}, {"CCA", 'P'}, {"CCG", 'P'},
        {"CAT", 'H'}, {"CAC", 'H'}, {"CAA", 'Q'}, {"CAG", 'Q'},
        {"CGT", 'R'}, {"CGC", 'R'}, {"CGA", 'R'}, {"CGG", 'R'},
        {"ATT", 'I'}, {"ATC", 'I'}, {"ATA", 'I'}, {"ATG", 'M'},
        {"ACT", 'T'}, {"ACC", 'T'}, {"ACA", 'T'}, {"ACG", 'T'},
        {"AAT", 'N'}, {"AAC", 'N'}, {"AAA", 'K'}, {"AAG", 'K'},
        {"AGT", 'S'}, {"AGC", 'S'}, {"AGA", 'R'}, {"AGG", 'R'},
        {"GTT", 'V'}, {"GTC", 'V'}, {"GTA", 'V'}, {"GTG", 'V'},
        {"GCT", 'A'}, {"GCC", 'A'}, {"GCA", 'A'}, {"GCG", 'A'},
        {"GAT", 'D'}, {"GAC", 'D'}, {"GAA", 'E'}, {"GAG", 'E'},
        {"GGT", 'G'}, {"GGC", 'G'}, {"GGA", 'G'}, {"GGG", 'G'},
    };
    return table.at(codon);
}

std::array<int, 64> make_codon_map() {
    std::array<int, 64> result{};
    for (int first = 0; first < 4; ++first) {
        for (int second = 0; second < 4; ++second) {
            for (int third = 0; third < 4; ++third) {
                std::string codon;
                codon.push_back(BASE_ORDER[first]);
                codon.push_back(BASE_ORDER[second]);
                codon.push_back(BASE_ORDER[third]);
                const char amino_acid = amino_acid_for_codon(codon);
                const int codon_index = first * 16 + second * 4 + third;
                if (amino_acid == '*') {
                    result[codon_index] = -1;
                } else {
                    const size_t amino_acid_index = AA_ORDER.find(amino_acid);
                    if (amino_acid_index == std::string::npos) {
                        throw std::runtime_error("Amino acid absent from AA_ORDER");
                    }
                    result[codon_index] = static_cast<int>(amino_acid_index);
                }
            }
        }
    }
    return result;
}

std::array<int, AA_COUNT> make_degeneracy(const std::array<int, 64>& codon_map) {
    std::array<int, AA_COUNT> degeneracy{};
    for (const int index : codon_map) if (index >= 0) ++degeneracy[index];
    const int total = std::accumulate(degeneracy.begin(), degeneracy.end(), 0);
    if (total != 61) throw std::runtime_error("Codon degeneracy does not sum to 61");
    return degeneracy;
}

TranslationCounts translated_counts(
    const std::array<std::string, 3>& sequences,
    bool all_forward_frames,
    const std::array<int, 64>& codon_map
) {
    TranslationCounts result;
    for (const auto& sequence : sequences) {
        const int frame_count = all_forward_frames ? 3 : 1;
        for (int frame = 0; frame < frame_count; ++frame) {
            for (size_t position = static_cast<size_t>(frame); position + 2 < sequence.size(); position += 3) {
                const int codon_index = base_index(sequence[position]) * 16
                    + base_index(sequence[position + 1]) * 4
                    + base_index(sequence[position + 2]);
                const int amino_acid_index = codon_map[codon_index];
                if (amino_acid_index < 0) {
                    ++result.stops;
                } else {
                    ++result.amino_acids[amino_acid_index];
                }
            }
        }
    }
    return result;
}

TranslationCounts protein_counts(const std::string& sequence) {
    TranslationCounts result;
    for (const char amino_acid : sequence) {
        const size_t index = AA_ORDER.find(amino_acid);
        if (index == std::string::npos) throw std::runtime_error("Unexpected protein residue");
        ++result.amino_acids[index];
    }
    return result;
}

uint32_t sense_total(const TranslationCounts& counts) {
    return std::accumulate(counts.amino_acids.begin(), counts.amino_acids.end(), uint32_t{0});
}

std::array<double, AA_COUNT> distribution(
    const TranslationCounts& counts,
    bool codon_normalized,
    const std::array<int, AA_COUNT>& degeneracy
) {
    std::array<double, AA_COUNT> result{};
    double denominator = 0.0;
    for (int index = 0; index < AA_COUNT; ++index) {
        result[index] = static_cast<double>(counts.amino_acids[index]);
        if (codon_normalized) result[index] /= static_cast<double>(degeneracy[index]);
        denominator += result[index];
    }
    if (denominator == 0.0) throw std::runtime_error("No sense amino acids in translation");
    for (double& value : result) value /= denominator;
    return result;
}

double js_divergence(
    const std::array<double, AA_COUNT>& first,
    const std::array<double, AA_COUNT>& second
) {
    double result = 0.0;
    for (int index = 0; index < AA_COUNT; ++index) {
        const double midpoint = 0.5 * (first[index] + second[index]);
        if (first[index] > 0.0) result += 0.5 * first[index] * std::log2(first[index] / midpoint);
        if (second[index] > 0.0) result += 0.5 * second[index] * std::log2(second[index] / midpoint);
    }
    return result;
}

std::array<std::string, 3> load_acro1_sequences(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot open ACRO1 input");
    std::string line;
    std::getline(input, line);
    const auto columns = header_index(line);
    std::map<std::string, std::string> sequences;
    while (std::getline(input, line)) {
        const auto fields = split_tab(line);
        sequences[fields.at(columns.at("allele"))] = fields.at(columns.at("species_balanced_consensus_nt"));
    }
    const std::array<std::string, 3> order = {
        "canonical_ACRO1", "internal_CCC", "boundary_terminal"
    };
    std::array<std::string, 3> result;
    for (size_t index = 0; index < order.size(); ++index) {
        result[index] = sequences.at(order[index]);
        if (result[index].find_first_not_of("ACGT") != std::string::npos) {
            throw std::runtime_error("Non-ACGT base in ACRO1 input");
        }
    }
    if (result[0].size() != 147 || result[1].size() != 144 || result[2].size() != 144) {
        throw std::runtime_error("Unexpected ACRO1 sequence lengths");
    }
    return result;
}

std::map<std::string, std::string> load_ctd_references(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot open POLR2A reference input");
    std::string line;
    std::getline(input, line);
    const auto columns = header_index(line);
    std::map<std::string, std::string> result;
    while (std::getline(input, line)) {
        const auto fields = split_tab(line);
        result[fields.at(columns.at("reference_id"))] = fields.at(columns.at("sequence"));
    }
    if (result.at("selected_segment_1728_1776").size() != 49) {
        throw std::runtime_error("Unexpected selected CTD segment length");
    }
    if (result.at("full_CTD_1593_1960").size() != 368) {
        throw std::runtime_error("Unexpected full CTD length");
    }
    const size_t offset = result.at("full_CTD_1593_1960").find(
        result.at("selected_segment_1728_1776")
    );
    if (offset != 135) throw std::runtime_error("Selected CTD segment is not at expected offset");
    return result;
}

std::array<double, 4> load_human_base_frequencies(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot open base-frequency input");
    std::string line;
    std::getline(input, line);
    const auto columns = header_index(line);
    std::map<char, double> frequencies;
    while (std::getline(input, line)) {
        const auto fields = split_tab(line);
        if (fields.at(columns.at("base")).size() != 1) continue;
        const char base = fields.at(columns.at("base"))[0];
        if (BASE_ORDER.find(base) == std::string::npos) continue;
        frequencies[base] = std::stod(fields.at(columns.at("fraction")));
    }
    std::array<double, 4> result{};
    double total = 0.0;
    for (int index = 0; index < 4; ++index) {
        result[index] = frequencies.at(BASE_ORDER[index]);
        total += result[index];
    }
    if (std::abs(total - 1.0) > 1e-12) throw std::runtime_error("Base frequencies do not sum to one");
    return result;
}

std::string human_iid_sequence(
    size_t length,
    const std::array<double, 4>& frequencies,
    std::mt19937_64& rng,
    std::uniform_real_distribution<double>& uniform
) {
    const double threshold_a = frequencies[0];
    const double threshold_c = threshold_a + frequencies[1];
    const double threshold_g = threshold_c + frequencies[2];
    std::string result(length, 'A');
    for (char& base : result) {
        const double draw = uniform(rng);
        if (draw < threshold_a) base = 'A';
        else if (draw < threshold_c) base = 'C';
        else if (draw < threshold_g) base = 'G';
        else base = 'T';
    }
    return result;
}

// Practical randomized Euler-trail dinucleotide shuffle. It exactly preserves
// sequence length, start/end bases, mononucleotide counts, and dinucleotide counts,
// but is not claimed to sample all Euler trails uniformly.
std::string dinucleotide_shuffle(const std::string& sequence, std::mt19937_64& rng) {
    if (sequence.size() < 2) return sequence;
    constexpr size_t MAX_LENGTH = 256;
    if (sequence.size() > MAX_LENGTH) throw std::runtime_error("Sequence exceeds shuffle buffer");
    std::array<std::array<uint8_t, MAX_LENGTH>, 4> adjacency{};
    std::array<int, 4> sizes{};
    for (size_t position = 0; position + 1 < sequence.size(); ++position) {
        const int source = base_index(sequence[position]);
        adjacency[source][sizes[source]++] = static_cast<uint8_t>(base_index(sequence[position + 1]));
    }
    for (int base = 0; base < 4; ++base) {
        std::shuffle(adjacency[base].begin(), adjacency[base].begin() + sizes[base], rng);
    }
    std::array<uint8_t, MAX_LENGTH> stack{};
    std::array<uint8_t, MAX_LENGTH> reverse_path{};
    int stack_size = 0;
    int path_size = 0;
    stack[stack_size++] = static_cast<uint8_t>(base_index(sequence.front()));
    while (stack_size > 0) {
        const uint8_t vertex = stack[stack_size - 1];
        if (sizes[vertex] > 0) {
            const uint8_t destination = adjacency[vertex][--sizes[vertex]];
            stack[stack_size++] = destination;
        } else {
            reverse_path[path_size++] = vertex;
            --stack_size;
        }
    }
    if (path_size != static_cast<int>(sequence.size())) {
        throw std::runtime_error("Dinucleotide Euler shuffle failed");
    }
    std::string result(sequence.size(), 'A');
    for (int index = 0; index < path_size; ++index) {
        result[index] = BASE_ORDER[reverse_path[path_size - 1 - index]];
    }
    return result;
}

double quantile(const std::vector<double>& sorted_values, double probability) {
    const double position = probability * static_cast<double>(sorted_values.size() - 1);
    const size_t lower = static_cast<size_t>(std::floor(position));
    const size_t upper = static_cast<size_t>(std::ceil(position));
    const double fraction = position - static_cast<double>(lower);
    return sorted_values[lower] * (1.0 - fraction) + sorted_values[upper] * fraction;
}

std::vector<TestAccumulator> initialize_tests(
    const TranslationCounts& observed_frame0,
    const TranslationCounts& observed_all_frames,
    const std::map<std::string, TranslationCounts>& ctd_counts,
    const std::array<int, AA_COUNT>& degeneracy,
    uint64_t replicates
) {
    std::vector<TestAccumulator> tests;
    for (const auto& frame_pair : std::vector<std::pair<std::string, TranslationCounts>>{
             {"candidate_frame0", observed_frame0},
             {"all_three_forward_frames", observed_all_frames},
         }) {
        for (const auto& ctd_pair : ctd_counts) {
            for (const bool codon_normalized : {false, true}) {
                const auto acro_distribution = distribution(frame_pair.second, codon_normalized, degeneracy);
                const auto ctd_distribution = distribution(ctd_pair.second, codon_normalized, degeneracy);
                TestAccumulator test;
                test.frame_scheme = frame_pair.first;
                test.ctd_reference = ctd_pair.first;
                test.normalization = codon_normalized ? "codon_degeneracy_normalized" : "raw_amino_acid_frequency";
                if (
                    frame_pair.first == "candidate_frame0"
                    && ctd_pair.first == "full_CTD_1593_1960"
                    && codon_normalized
                ) {
                    test.analysis_role = "candidate_ORF_primary";
                } else if (
                    frame_pair.first == "all_three_forward_frames"
                    && ctd_pair.first == "full_CTD_1593_1960"
                    && codon_normalized
                ) {
                    test.analysis_role = "strict_RAP_method_replication";
                } else {
                    test.analysis_role = "sensitivity_analysis";
                }
                test.observed_jsd = js_divergence(acro_distribution, ctd_distribution);
                test.observed_sense_residues = sense_total(frame_pair.second);
                test.observed_stop_codons = frame_pair.second.stops;
                test.values.reserve(static_cast<size_t>(replicates));
                tests.push_back(std::move(test));
            }
        }
    }
    return tests;
}

void add_simulated_values(
    std::vector<TestAccumulator>& tests,
    const TranslationCounts& simulated_frame0,
    const TranslationCounts& simulated_all_frames,
    const std::map<std::string, TranslationCounts>& ctd_counts,
    const std::array<int, AA_COUNT>& degeneracy
) {
    size_t test_index = 0;
    for (const auto& frame_pair : std::vector<std::pair<std::string, TranslationCounts>>{
             {"candidate_frame0", simulated_frame0},
             {"all_three_forward_frames", simulated_all_frames},
         }) {
        for (const auto& ctd_pair : ctd_counts) {
            for (const bool codon_normalized : {false, true}) {
                const auto acro_distribution = distribution(frame_pair.second, codon_normalized, degeneracy);
                const auto ctd_distribution = distribution(ctd_pair.second, codon_normalized, degeneracy);
                tests[test_index].values.push_back(js_divergence(acro_distribution, ctd_distribution));
                ++test_index;
            }
        }
    }
}

void write_composition_rows(
    std::ofstream& composition_output,
    const std::string& source_id,
    const std::string& source_type,
    const std::string& frame_scheme,
    const TranslationCounts& counts,
    const std::array<int, AA_COUNT>& degeneracy
) {
    const auto raw = distribution(counts, false, degeneracy);
    const auto normalized = distribution(counts, true, degeneracy);
    const uint32_t total = sense_total(counts);
    for (int index = 0; index < AA_COUNT; ++index) {
        composition_output << source_id << '\t' << source_type << '\t' << frame_scheme << '\t'
            << AA_ORDER[index] << '\t' << counts.amino_acids[index] << '\t' << total << '\t'
            << counts.stops << '\t' << raw[index] << '\t' << degeneracy[index] << '\t'
            << normalized[index] << '\n';
    }
}

void write_observed_summary_row(
    std::ofstream& output,
    const std::string& source_id,
    const std::string& source_type,
    const std::string& frame_scheme,
    const TranslationCounts& counts
) {
    const auto count_for = [&](char amino_acid) {
        return counts.amino_acids[AA_ORDER.find(amino_acid)];
    };
    const uint32_t total = sense_total(counts);
    const uint32_t s = count_for('S');
    const uint32_t p = count_for('P');
    const uint32_t t = count_for('T');
    const uint32_t y = count_for('Y');
    output << source_id << '\t' << source_type << '\t' << frame_scheme << '\t' << total << '\t'
        << counts.stops << '\t' << s << '\t' << p << '\t' << t << '\t' << y << '\t'
        << (s + p + t) << '\t' << (s + p + t + y) << '\t'
        << static_cast<double>(s + p + t) / total << '\t'
        << static_cast<double>(s + p + t + y) / total << '\n';
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 10) {
        std::cerr << "usage: jsd_null_simulator ACRO_TSV POLR2A_TSV BASE_TSV REPLICATES SEED SUMMARY_OUT HIST_OUT COMPOSITION_OUT OBSERVED_OUT\n";
        return 1;
    }
    try {
        const std::string acro_path = argv[1];
        const std::string polr2a_path = argv[2];
        const std::string base_path = argv[3];
        const uint64_t replicates = std::stoull(argv[4]);
        const uint64_t seed = std::stoull(argv[5]);
        const std::string summary_path = argv[6];
        const std::string histogram_path = argv[7];
        const std::string composition_path = argv[8];
        const std::string observed_path = argv[9];

        const auto codon_map = make_codon_map();
        const auto degeneracy = make_degeneracy(codon_map);
        const auto observed_sequences = load_acro1_sequences(acro_path);
        const auto ctd_sequences = load_ctd_references(polr2a_path);
        const auto human_base_frequencies = load_human_base_frequencies(base_path);

        std::map<std::string, TranslationCounts> ctd_counts;
        for (const auto& pair : ctd_sequences) ctd_counts[pair.first] = protein_counts(pair.second);
        const TranslationCounts observed_frame0 = translated_counts(observed_sequences, false, codon_map);
        const TranslationCounts observed_all_frames = translated_counts(observed_sequences, true, codon_map);

        std::ofstream summary_output(summary_path);
        std::ofstream histogram_output(histogram_path);
        std::ofstream composition_output(composition_path);
        std::ofstream observed_output(observed_path);
        if (!summary_output || !histogram_output || !composition_output || !observed_output) {
            throw std::runtime_error("Cannot open one or more output files");
        }
        summary_output << std::setprecision(17);
        histogram_output << std::setprecision(17);
        composition_output << std::setprecision(17);
        observed_output << std::setprecision(17);

        summary_output
            << "analysis_role\tnull_model\tnull_model_definition\tframe_scheme\tctd_reference\tnormalization"
            << "\tobserved_jsd_base2\tobserved_sense_residues\tobserved_stop_codons\tnull_replicates"
            << "\tnull_jsd_less_equal_observed\tempirical_p_plus1\tzero_lower_tail_event_95_upper_bound"
            << "\tnull_mean\tnull_sd\tnull_min\tnull_q025\tnull_median\tnull_q975\tnull_max"
            << "\tstandardized_difference\tobserved_lower_tail_percentile\trandom_seed\n";
        histogram_output
            << "null_model\tframe_scheme\tctd_reference\tnormalization\tbin_start\tbin_end\tcount\n";
        composition_output
            << "source_id\tsource_type\tframe_scheme\tamino_acid\tcount\tsense_total\tstop_count"
            << "\traw_fraction\tcodon_degeneracy\tcodon_normalized_fraction\n";
        observed_output
            << "source_id\tsource_type\tframe_scheme\tsense_residues\tstop_codons\tS_count\tP_count"
            << "\tT_count\tY_count\tSPT_count\tSPTY_count\tSPT_fraction\tSPTY_fraction\n";

        write_composition_rows(
            composition_output, "ACRO1_species_balanced_pooled_frame0", "ACRO1", "candidate_frame0",
            observed_frame0, degeneracy
        );
        write_composition_rows(
            composition_output, "ACRO1_species_balanced_pooled_all_forward_frames", "ACRO1",
            "all_three_forward_frames", observed_all_frames, degeneracy
        );
        write_observed_summary_row(
            observed_output, "ACRO1_species_balanced_pooled_frame0", "ACRO1", "candidate_frame0",
            observed_frame0
        );
        write_observed_summary_row(
            observed_output, "ACRO1_species_balanced_pooled_all_forward_frames", "ACRO1",
            "all_three_forward_frames", observed_all_frames
        );
        for (size_t index = 0; index < observed_sequences.size(); ++index) {
            std::array<std::string, 3> single = {observed_sequences[index], "", ""};
            const TranslationCounts counts = translated_counts(single, false, codon_map);
            const std::array<std::string, 3> names = {
                "ACRO1_canonical_frame0", "ACRO1_internal_CCC_frame0", "ACRO1_boundary_terminal_frame0"
            };
            write_observed_summary_row(observed_output, names[index], "ACRO1_allele", "candidate_frame0", counts);
        }
        for (const auto& pair : ctd_counts) {
            write_composition_rows(composition_output, pair.first, "POLR2A_CTD", "not_applicable", pair.second, degeneracy);
            write_observed_summary_row(observed_output, pair.first, "POLR2A_CTD", "not_applicable", pair.second);
        }

        const std::vector<std::pair<std::string, std::string>> null_models = {
            {
                "human_genome_iid",
                "Independent nucleotides sampled from whole-HG002 A/C/G/T frequencies; sequence lengths matched."
            },
            {
                "exact_mononucleotide_shuffle",
                "Each ACRO1 consensus independently shuffled; exact sequence-specific nucleotide counts preserved."
            },
            {
                "dinucleotide_shuffle",
                "Randomized Euler-trail shuffle; exact sequence-specific mono- and dinucleotide counts preserved."
            },
        };

        std::uniform_real_distribution<double> uniform(0.0, 1.0);
        for (size_t null_index = 0; null_index < null_models.size(); ++null_index) {
            const std::string& null_model = null_models[null_index].first;
            const std::string& null_definition = null_models[null_index].second;
            const uint64_t model_seed = seed + static_cast<uint64_t>(null_index) * 1000003ULL;
            std::mt19937_64 rng(model_seed);
            auto tests = initialize_tests(
                observed_frame0, observed_all_frames, ctd_counts, degeneracy, replicates
            );

            for (uint64_t replicate = 0; replicate < replicates; ++replicate) {
                std::array<std::string, 3> simulated_sequences;
                for (size_t sequence_index = 0; sequence_index < observed_sequences.size(); ++sequence_index) {
                    if (null_model == "human_genome_iid") {
                        simulated_sequences[sequence_index] = human_iid_sequence(
                            observed_sequences[sequence_index].size(), human_base_frequencies, rng, uniform
                        );
                    } else if (null_model == "exact_mononucleotide_shuffle") {
                        simulated_sequences[sequence_index] = observed_sequences[sequence_index];
                        std::shuffle(
                            simulated_sequences[sequence_index].begin(),
                            simulated_sequences[sequence_index].end(),
                            rng
                        );
                    } else {
                        simulated_sequences[sequence_index] = dinucleotide_shuffle(
                            observed_sequences[sequence_index], rng
                        );
                    }
                }
                const TranslationCounts simulated_frame0 = translated_counts(
                    simulated_sequences, false, codon_map
                );
                const TranslationCounts simulated_all_frames = translated_counts(
                    simulated_sequences, true, codon_map
                );
                add_simulated_values(
                    tests, simulated_frame0, simulated_all_frames, ctd_counts, degeneracy
                );
            }

            for (auto& test : tests) {
                const double sum = std::accumulate(test.values.begin(), test.values.end(), 0.0);
                const double mean = sum / static_cast<double>(test.values.size());
                double squared_sum = 0.0;
                for (const double value : test.values) squared_sum += (value - mean) * (value - mean);
                const double standard_deviation = std::sqrt(
                    squared_sum / static_cast<double>(test.values.size())
                );
                std::sort(test.values.begin(), test.values.end());
                const uint64_t less_equal = static_cast<uint64_t>(
                    std::upper_bound(test.values.begin(), test.values.end(), test.observed_jsd)
                    - test.values.begin()
                );
                const double empirical_p = static_cast<double>(less_equal + 1)
                    / static_cast<double>(replicates + 1);
                const double zero_upper = less_equal == 0
                    ? 1.0 - std::pow(0.05, 1.0 / static_cast<double>(replicates))
                    : std::numeric_limits<double>::quiet_NaN();
                const double standardized_difference = (test.observed_jsd - mean) / standard_deviation;
                summary_output
                    << test.analysis_role << '\t' << null_model << '\t' << null_definition << '\t'
                    << test.frame_scheme << '\t' << test.ctd_reference << '\t' << test.normalization << '\t'
                    << test.observed_jsd << '\t' << test.observed_sense_residues << '\t'
                    << test.observed_stop_codons << '\t' << replicates << '\t' << less_equal << '\t'
                    << empirical_p << '\t';
                if (std::isnan(zero_upper)) summary_output << "";
                else summary_output << zero_upper;
                summary_output
                    << '\t' << mean << '\t' << standard_deviation << '\t' << test.values.front()
                    << '\t' << quantile(test.values, 0.025) << '\t' << quantile(test.values, 0.5)
                    << '\t' << quantile(test.values, 0.975) << '\t' << test.values.back()
                    << '\t' << standardized_difference << '\t'
                    << static_cast<double>(less_equal) / static_cast<double>(replicates)
                    << '\t' << model_seed << '\n';

                std::array<uint64_t, HISTOGRAM_BINS> histogram{};
                for (const double value : test.values) {
                    int bin = static_cast<int>(std::floor(value * HISTOGRAM_BINS));
                    if (bin < 0) bin = 0;
                    if (bin >= HISTOGRAM_BINS) bin = HISTOGRAM_BINS - 1;
                    ++histogram[bin];
                }
                for (int bin = 0; bin < HISTOGRAM_BINS; ++bin) {
                    histogram_output << null_model << '\t' << test.frame_scheme << '\t'
                        << test.ctd_reference << '\t' << test.normalization << '\t'
                        << static_cast<double>(bin) / HISTOGRAM_BINS << '\t'
                        << static_cast<double>(bin + 1) / HISTOGRAM_BINS << '\t'
                        << histogram[bin] << '\n';
                }
            }
            std::cerr << "completed " << null_model << " with " << replicates << " replicates\n";
        }
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 2;
    }
    return 0;
}
