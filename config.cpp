#include "config.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <string>

qual_t Config::Q_MIN;   /* Minimum quality threshold. */
count_t Config::D_MIN;   /* Minimum number of high-quality extensions needed. */

k_t Config::K;

size_t Config::MIN_CONTIG_LEN;
size_t Config::FASTA_TEXTWIDTH;

void Config::load_config (const std::string &config_path) {

    boost::property_tree::ptree pTree;

    try {
        read_info(config_path, pTree);
    }
    catch (boost::property_tree::info_parser_error e) {
        std::cout << "Error reading config file." << std::endl;
    }

    /* Check if defaults need to be overwritten */
    Config::Q_MIN = pTree.get("defaults.Q_MIN", 19); /* Minimum quality threshold */
    Config::D_MIN = pTree.get("defaults.D_MIN", 10); /* Minimum number of high-quality extensions needed. */
    Config::K = pTree.get("defaults.K", 41);
    Config::MIN_CONTIG_LEN = pTree.get("defaults.MIN_CONTIG_LEN", 2*K);
    Config::FASTA_TEXTWIDTH = pTree.get("defaults.FASTA_TEXTWIDTH", 50);

}
