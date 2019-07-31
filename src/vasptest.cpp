#include <cassert>

#include "al_VASPReader.hpp"

int main(int argc, char *argv[])
{
    VASPReader reader;
    reader.loadFile("/2809share/casm_viewer/bin/cached_output/_alloshare_vdv group_proj1_monte_heating_template_POSCAR_a_monte_heatingchempot_-1_050_conditions_29_");


    VASPReader reader2;
    reader2.loadFile("/2809share/casm_viewer/bin/cached_output/_alloshare_vdv group_proj1_monte_heating_template_POSCAR_a_monte_heatingchempot_-1_050_conditions_29_");

    auto &elems1 = reader.getAllPositions();
    auto &elems2 = reader2.getAllPositions();

    assert(elems1.size() == elems2.size());

    for (auto &atoms: elems1) {
        auto atoms2 = reader2.getElementPositions(atoms.first);
        assert(atoms2.size() == atoms.second.size());
    }

    return 0;

}
