# -*- mode: cmake -*-

#
#  ATS
#    internal energy constitutive relations
#

add_library(energy_relations_enthalpy
  enthalpy_evaluator.cc
             )

install(TARGETS energy_relations_enthalpy DESTINATION lib)
