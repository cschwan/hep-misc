#ifdef RECOLA
#define GNU_NAME_MANGLING
#include <recola.h>
#endif

#include <array>
#include <cassert>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

enum conv_t {
    blha = 0,
    coli = 1
};

#ifdef OPENLOOPS
extern "C"
{

    void ol_setparameter_int(const char* param, int val);
    void ol_setparameter_double(const char* param, double val);
    int ol_register_process(const char* process, int amptype);
    void ol_start();
    void ol_finish();
    void ol_evaluate_tree(int id, double* pp, double* m2tree);
    void ol_evaluate_cc(int id, double* pp, double* m2tree, double* m2cc, double* m2ew);
    void ol_evaluate_loop(int id, double* pp, double* m2tree, double* m2loop, double* acc);
}

int openloops_1(int alphas, conv_t conv) {
    // TODO: additional parameters not set for Recola
    ol_setparameter_int("ew_scheme", 1);
    ol_setparameter_int("ew_renorm_scheme", 1);
    ol_setparameter_int("ew_renorm", 1);
    ol_setparameter_int("verbose", 2);

    ol_setparameter_double("mass(23)", 91.153480619183);
    ol_setparameter_double("width(23)", 2.4942663787728);
    ol_setparameter_double("mass(24)", 80.351984549666);
    ol_setparameter_double("width(24)", 2.0837993977753);
    ol_setparameter_double("mass(25)", 125.0);
    ol_setparameter_double("width(25)", 4.07e-3);
    ol_setparameter_double("mass(6)", 173.0);
    ol_setparameter_double("width(6)", 0.0);
    ol_setparameter_double("gmu", 1.1663787e-5);
    ol_setparameter_int("complex_mass_scheme", 2);

    ol_setparameter_double("alphas", 1.0);
    ol_setparameter_double("mureg", 100.0);
    ol_setparameter_double("muren", 100.0);

    //ol_setparameter_int("write_parameters", 1);
    //ol_setparameter_int("verbose", 1);

    switch (conv)
    {
    case blha:
        ol_setparameter_int("polenorm", 0);
        break;

    case coli:
        ol_setparameter_int("polenorm", 1);
        break;

    default:
        assert(false);
    }

    // color-correlated matrix elements are in the library containing the virtuals
    ol_setparameter_int("loop_order_qcd", alphas);
    int id = ol_register_process("2 2 -> -11 12 -13 14 1 1", 11);

    //ol_setparameter_int("order_qcd", alphas);
    //int id = ol_register_process("2 2 -> -11 12 -13 14 1 1", 1);

    ol_start();

    return id;
}

double openloops_2(int id, double mom[]) {
    // double result;
    // ol_evaluate_tree(id, mom, &result);
    // return result;

    double m2tree;
    //std::vector<double> m2cc(8 * (8 - 1) / 2);
    //double m2ew;
    //ol_evaluate_cc(id, mom, &m2tree, m2cc.data(), &m2ew);
    //double result = m2cc.at(i + j * (j - 1) / 2);
    double m2loop[3];
    double acc;
    ol_evaluate_loop(id, mom, &m2tree, m2loop, &acc);

    //std::cout << "m[0]: " << m2loop[0] << '\n';
    //std::cout << "m[1]: " << m2loop[1] << '\n';
    //std::cout << "m[2]: " << m2loop[2] << '\n';
    //std::cout << "mtree:" << m2tree << '\n';

    //result = m2loop[0];
    //double pi = std::acos(-1.0);
    //result += pi * pi / 6.0 * m2loop[2];

    return m2loop[0];
}
#endif

#ifdef RECOLA
int recola_1(int alphas, conv_t conv) {
    Recola::set_complex_mass_scheme_rcl();
    Recola::set_pole_mass_z_rcl(91.153480619183, 2.4942663787728);
    Recola::set_pole_mass_w_rcl(80.351984549666, 2.0837993977753);
    Recola::set_pole_mass_h_rcl(125.0, 4.07e-3);
    Recola::set_pole_mass_top_rcl(173.0, 0.0);
    Recola::use_gfermi_scheme_real_and_set_gfermi_rcl(1.1663787e-5);
    //Recola::use_gfermi_scheme_and_set_gfermi_rcl(1.1663787e-5);
    Recola::set_print_level_squared_amplitude_rcl(1);

    switch (conv)
    {
    case blha:
        {
            double pi = std::acos(-1.0);
            Recola::set_delta_ir_rcl(0.0, pi * pi / 6.0);
        }
        break;

    case coli:
        break;

    default:
        assert(false);
    }

    // this isn't Q=100 GeV, but we don't vary the scale so it shouldn't matter
    Recola::set_alphas_rcl(1.0, 100.0, 5);

    int id = 1;

    Recola::define_process_rcl(id, "u u -> e+ nu_e mu+ nu_mu d d", "NLO");
    // Recola::define_process_rcl(1, "u u -> e+ nu_e mu+ nu_mu d d", "LO");
    Recola::unselect_all_gs_powers_LoopAmpl_rcl(id);
    Recola::unselect_all_gs_powers_BornAmpl_rcl(id);

    // Recola::select_gs_power_BornAmpl_rcl(1, alphas);

    std::array<int, 2> born_gs = { 0, 2 };
    std::array<int, 3> loop_gs = { 0, 2, 4 };

    for (int born : born_gs) {
        for (int loop : loop_gs) {
            if ((born + loop) == 2 * alphas) {
                Recola::select_gs_power_LoopAmpl_rcl(id, loop);
                Recola::select_gs_power_BornAmpl_rcl(id, born);
            }
        }
    }

    Recola::generate_processes_rcl();

    return 1;
}

double recola_2(int id, double alphas, double mom[][4]) {
    // double result = 0.0;
    // Recola::set_print_level_squared_amplitude_rcl(1);
    // Recola::compute_process_rcl(1, mom, "LO");
    // Recola::get_squared_amplitude_rcl(1, alphas, "LO", result);

    // return result;

    //double result = 0.0;
    //Recola::set_print_level_correlations_rcl(1);
    //Recola::compute_colour_correlation_rcl(1, mom, i + 1, j + 1);
    //Recola::get_colour_correlation_rcl(1, alphas - 1, i + 1, j + 1, result);
    //result *= 4.0 / 3.0;
    double result[2] = {};
    Recola::compute_process_rcl(1, mom, "NLO", result);
    Recola::get_squared_amplitude_rcl(1, alphas, "NLO", result[1]);

    //std::cout << "a(0): " << result[0] << '\n';
    //std::cout << "a(1): " << result[1] << '\n';
    //std::cout << "a   : " << result[2] << '\n';

    return result[1];
}
#endif

int main(int argc, char* argv[]) {
    int alphas = std::stoi(argv[1]);
    conv_t conv = static_cast <conv_t> (std::stoi(argv[2]));

#ifdef RECOLA
    int rid = recola_1(alphas, conv);
#endif

#ifdef OPENLOOPS
    int oid = openloops_1(alphas, conv);
#endif

    std::ifstream input(argv[3]);
    std::string line;

    std::cout << "O(as^" << alphas << ") ";
    switch (conv) {
    case blha:
        std::cout << "(BLHA): \n";
        break;

    case coli:
        std::cout << "(COLI): \n";
        break;

    default:
        assert(false);
    }

    while (std::getline(input, line))
    {
        std::istringstream in(line);
        std::vector<double> x(32);

        for (std::size_t i = 0; i != x.size(); ++i)
        {
            if (!(in >> x[i]))
            {
                std::cerr << "Error while reading input\n";
                return 1;
            }
        }

#ifdef RECOLA
        double rmom[8][4] = {
            { x[ 0], x[ 1], x[ 2], x[ 3] },
            { x[ 4], x[ 5], x[ 6], x[ 7] },
            { x[ 8], x[ 9], x[10], x[11] },
            { x[12], x[13], x[14], x[15] },
            { x[16], x[17], x[18], x[19] },
            { x[20], x[21], x[22], x[23] },
            { x[24], x[25], x[26], x[27] },
            { x[28], x[29], x[30], x[31] },
        };
        double recola = recola_2(rid, alphas, rmom);
        std::cout << std::setprecision(13) << recola << std::endl;
#endif

#ifdef OPENLOOPS
        double mom[40] = {
            x[ 0], x[ 1], x[ 2], x[ 3], 0.0,
            x[ 4], x[ 5], x[ 6], x[ 7], 0.0,
            x[ 8], x[ 9], x[10], x[11], 0.0,
            x[12], x[13], x[14], x[15], 0.0,
            x[16], x[17], x[18], x[19], 0.0,
            x[20], x[21], x[22], x[23], 0.0,
            x[24], x[25], x[26], x[27], 0.0,
            x[28], x[29], x[30], x[31], 0.0,
        };

        double openloops = openloops_2(oid, mom);
        std::cout << std::setprecision(13) << openloops << std::endl;
#endif
    }

#ifdef RECOLA
#endif

#ifdef OPENLOOPS
    ol_finish();
#endif
}
