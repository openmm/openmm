from typing import Union


class Units:
    HARTREE_TO_KJMOL = 2625.4996
    KJMOL_TO_HARTREE = 1.0 / HARTREE_TO_KJMOL

    BOHR_TO_NM = 0.052917721
    NM_TO_BOHR = 1.0 / BOHR_TO_NM

    AU_TIME_TO_PS = 2.4188843265e-5
    PS_TO_AU_TIME = 1.0 / AU_TIME_TO_PS

    AMU_TO_AU_MASS = 1822.888
    AU_MASS_TO_AMU = 1.0 / AMU_TO_AU_MASS

    HARTREE_TO_CM1 = 219474.63
    CM1_TO_HARTREE = 1.0 / HARTREE_TO_CM1

    KB_KJMOL_PER_K = 0.008314462618
    KB_HARTREE_PER_K = 3.16681153e-6

    @staticmethod
    def cm1_to_au(freq_cm1: float) -> float:
        return freq_cm1 * Units.CM1_TO_HARTREE

    @staticmethod
    def au_to_cm1(freq_au: float) -> float:
        return freq_au * Units.HARTREE_TO_CM1

    @staticmethod
    def kelvin_to_kT_kjmol(T_kelvin: float) -> float:
        return T_kelvin * Units.KB_KJMOL_PER_K

    @staticmethod
    def kT_kjmol_to_kelvin(kT_kjmol: float) -> float:
        return kT_kjmol / Units.KB_KJMOL_PER_K

    @staticmethod
    def kelvin_to_kT_hartree(T_kelvin: float) -> float:
        return T_kelvin * Units.KB_HARTREE_PER_K

    @staticmethod
    def kT_hartree_to_kelvin(kT_hartree: float) -> float:
        return kT_hartree / Units.KB_HARTREE_PER_K

    @staticmethod
    def kjmol_to_hartree(energy_kjmol: float) -> float:
        return energy_kjmol * Units.KJMOL_TO_HARTREE

    @staticmethod
    def hartree_to_kjmol(energy_hartree: float) -> float:
        return energy_hartree * Units.HARTREE_TO_KJMOL
