import pandas as pd
import os
import re
from collections import Counter
import numpy as np
from numba import jit

class Essentials:
    IUPAC = {
        "A" : ["A"], 
        "C" : ["C"], 
        "G" : ["G"],        
        "T" : ["T"],        
        "U" : ["U"],        
        "R" : ["G", "A"], 
        "Y" : ["T", "C"], 
        "K" : ["G", "T"], 
        "M" : ["A", "C"], 
        "S" : ["G", "C"], 
        "W" : ["A", "T"], 
        "B" : ["C", "G", "T"], 
        "D" : ["A", "G", "T"], 
        "H" : ["A", "C", "T"], 
        "V" : ["A", "C", "G"], 
        "N" : ["A", "C", "G", "T"],
        "-" : ["-"]
    }
    COMPLEMENTARY = {
        "A" : "T",
        "C" : "G",
        "G" : "C",
        "T" : "A",
        "R" : "Y",
        "Y" : "R",
        "K" : "M",
        "M" : "K",
        "S" : "S",
        "W" : "W",
        "B" : "V",
        "D" : "H",
        "H" : "D",
        "V" : "B",
        "N" : "N"   
    }
    def Antisense(sequence):
        nuclist = list(sequence)
        comp_nuclist = nuclist.copy()
        for i in range(len(nuclist)):
            comp_nuclist[i] = Essentials.COMPLEMENTARY[nuclist[i]]
        return ("").join(comp_nuclist)

    def GCcontent(sequence):
        gc = sequence.count("G") + sequence.count("C")
        gc_content = gc / len(sequence)
        return gc_content * 100
        
    def Tm(seq,
           GC,
           Na = 50,
           K = 0,
           Tris = 0,
           Mg = 0,
           dNTPs = 0,
           shift = 0,
           nn_table = None,
           tmm_table = None,
           imm_table = None,
           de_table = None,
           dnac1 = 25,
           dnac2 = 25,
           salt_correction = False):
        ############ Melting temperature calculations were cloned from BioPython library
        # https://github.com/bioinf-mcb/biopython/blob/master/Bio/SeqUtils/MeltingTemp.py

        # Thermodynamic lookup tables (dictionaries):
        # Enthalpy (dH) and entropy (dS) values for nearest neighbors and initiation
        # process. Calculation of duplex initiation is quite different in several
        # papers; to allow for a general calculation, all different initiation
        # parameters are included in all tables and non-applicable parameters are set
        # to zero.
        # The key is either an initiation type (e.g., 'init_A/T') or a nearest neighbor
        # duplex sequence (e.g., GT/CA, to read 5'GT3'-3'CA5'). The values are tuples
        # of dH (kcal/mol), dS (cal/mol K).

        # Turn black code style off
        # fmt: off

        # DNA/DNA
        # Breslauer et al. (1986), Proc Natl Acad Sci USA 83: 3746-3750
        ############ Melting temperature calculations were cloned from BioPython library
        # https://github.com/bioinf-mcb/biopython/blob/master/Bio/SeqUtils/MeltingTemp.py

        # Thermodynamic lookup tables (dictionaries):
        # Enthalpy (dH) and entropy (dS) values for nearest neighbors and initiation
        # process. Calculation of duplex initiation is quite different in several
        # papers; to allow for a general calculation, all different initiation
        # parameters are included in all tables and non-applicable parameters are set
        # to zero.
        # The key is either an initiation type (e.g., 'init_A/T') or a nearest neighbor
        # duplex sequence (e.g., GT/CA, to read 5'GT3'-3'CA5'). The values are tuples
        # of dH (kcal/mol), dS (cal/mol K).

        # Turn black code style off
        # fmt: off

        # DNA/DNA
        # Breslauer et al. (1986), Proc Natl Acad Sci USA 83: 3746-3750
        DNA_NN1 = {
            "init": (0, 0), "init_A/T": (0, 0), "init_G/C": (0, 0),
            "init_oneG/C": (0, -16.8), "init_allA/T": (0, -20.1), "init_5T/A": (0, 0),
            "sym": (0, -1.3),
            "AA/TT": (-9.1, -24.0), "AT/TA": (-8.6, -23.9), "TA/AT": (-6.0, -16.9),
            "CA/GT": (-5.8, -12.9), "GT/CA": (-6.5, -17.3), "CT/GA": (-7.8, -20.8),
            "GA/CT": (-5.6, -13.5), "CG/GC": (-11.9, -27.8), "GC/CG": (-11.1, -26.7),
            "GG/CC": (-11.0, -26.6)}

        # Sugimoto et al. (1996), Nuc Acids Res 24 : 4501-4505
        DNA_NN2 = {
            "init": (0.6, -9.0), "init_A/T": (0, 0), "init_G/C": (0, 0),
            "init_oneG/C": (0, 0), "init_allA/T": (0, 0), "init_5T/A": (0, 0),
            "sym": (0, -1.4),
            "AA/TT": (-8.0, -21.9), "AT/TA": (-5.6, -15.2), "TA/AT": (-6.6, -18.4),
            "CA/GT": (-8.2, -21.0), "GT/CA": (-9.4, -25.5), "CT/GA": (-6.6, -16.4),
            "GA/CT": (-8.8, -23.5), "CG/GC": (-11.8, -29.0), "GC/CG": (-10.5, -26.4),
            "GG/CC": (-10.9, -28.4)}

        # Allawi and SantaLucia (1997), Biochemistry 36: 10581-10594
        DNA_NN3 = {
            "init": (0, 0), "init_A/T": (2.3, 4.1), "init_G/C": (0.1, -2.8),
            "init_oneG/C": (0, 0), "init_allA/T": (0, 0), "init_5T/A": (0, 0),
            "sym": (0, -1.4),
            "AA/TT": (-7.9, -22.2), "AT/TA": (-7.2, -20.4), "TA/AT": (-7.2, -21.3),
            "CA/GT": (-8.5, -22.7), "GT/CA": (-8.4, -22.4), "CT/GA": (-7.8, -21.0),
            "GA/CT": (-8.2, -22.2), "CG/GC": (-10.6, -27.2), "GC/CG": (-9.8, -24.4),
            "GG/CC": (-8.0, -19.9)}

        # SantaLucia & Hicks (2004), Annu. Rev. Biophys. Biomol. Struct 33: 415-440
        DNA_NN4 = {
            "init": (0.2, -5.7), "init_A/T": (2.2, 6.9), "init_G/C": (0, 0),
            "init_oneG/C": (0, 0), "init_allA/T": (0, 0), "init_5T/A": (0, 0),
            "sym": (0, -1.4),
            "AA/TT": (-7.6, -21.3), "AT/TA": (-7.2, -20.4), "TA/AT": (-7.2, -20.4),
            "CA/GT": (-8.5, -22.7), "GT/CA": (-8.4, -22.4), "CT/GA": (-7.8, -21.0),
            "GA/CT": (-8.2, -22.2), "CG/GC": (-10.6, -27.2), "GC/CG": (-9.8, -24.4),
            "GG/CC": (-8.0, -19.0)}

        # RNA/RNA
        # Freier et al. (1986), Proc Natl Acad Sci USA 83: 9373-9377
        RNA_NN1 = {
            "init": (0, -10.8), "init_A/T": (0, 0), "init_G/C": (0, 0),
            "init_oneG/C": (0, 0), "init_allA/T": (0, 0), "init_5T/A": (0, 0),
            "sym": (0, -1.4),
            "AA/TT": (-6.6, -18.4), "AT/TA": (-5.7, -15.5), "TA/AT": (-8.1, -22.6),
            "CA/GT": (-10.5, -27.8), "GT/CA": (-10.2, -26.2), "CT/GA": (-7.6, -19.2),
            "GA/CT": (-13.3, -35.5), "CG/GC": (-8.0, -19.4), "GC/CG": (-14.2, -34.9),
            "GG/CC": (-12.2, -29.7)}

        # Xia et al (1998), Biochemistry 37: 14719-14735
        RNA_NN2 = {
            "init": (3.61, -1.5), "init_A/T": (3.72, 10.5), "init_G/C": (0, 0),
            "init_oneG/C": (0, 0), "init_allA/T": (0, 0), "init_5T/A": (0, 0),
            "sym": (0, -1.4),
            "AA/TT": (-6.82, -19.0), "AT/TA": (-9.38, -26.7), "TA/AT": (-7.69, -20.5),
            "CA/GT": (-10.44, -26.9), "GT/CA": (-11.40, -29.5),
            "CT/GA": (-10.48, -27.1), "GA/CT": (-12.44, -32.5),
            "CG/GC": (-10.64, -26.7), "GC/CG": (-14.88, -36.9),
            "GG/CC": (-13.39, -32.7)}

        # Chen et al. (2012), Biochemistry 51: 3508-3522
        RNA_NN3 = {
            "init": (6.40, 6.99), "init_A/T": (3.85, 11.04), "init_G/C": (0, 0),
            "init_oneG/C": (0, 0), "init_allA/T": (0, 0), "init_5T/A": (0, 0),
            "sym": (0, -1.4),
            "AA/TT": (-7.09, -19.8), "AT/TA": (-9.11, -25.8), "TA/AT": (-8.50, -22.9),
            "CA/GT": (-11.03, -28.8), "GT/CA": (-11.98, -31.3),
            "CT/GA": (-10.90, -28.5), "GA/CT": (-13.21, -34.9),
            "CG/GC": (-10.88, -27.4), "GC/CG": (-16.04, -40.6),
            "GG/CC": (-14.18, -35.0), "GT/TG": (-13.83, -46.9),
            "GG/TT": (-17.82, -56.7), "AG/TT": (-3.96, -11.6),
            "TG/AT": (-0.96, -1.8), "TT/AG": (-10.38, -31.8), "TG/GT": (-12.64, -38.9),
            "AT/TG": (-7.39, -21.0), "CG/GT": (-5.56, -13.9), "CT/GG": (-9.44, -24.7),
            "GG/CT": (-7.03, -16.8), "GT/CG": (-11.09, -28.8)}

        # RNA/DNA
        # Sugimoto et al. (1995), Biochemistry 34: 11211-11216
        R_DNA_NN1 = {
            "init": (1.9, -3.9), "init_A/T": (0, 0), "init_G/C": (0, 0),
            "init_oneG/C": (0, 0), "init_allA/T": (0, 0), "init_5T/A": (0, 0),
            "sym": (0, 0),
            "AA/TT": (-11.5, -36.4), "AC/TG": (-7.8, -21.6), "AG/TC": (-7.0, -19.7),
            "AT/TA": (-8.3, -23.9), "CA/GT": (-10.4, -28.4), "CC/GG": (-12.8, -31.9),
            "CG/GC": (-16.3, -47.1), "CT/GA": (-9.1, -23.5), "GA/CT": (-8.6, -22.9),
            "GC/CG": (-8.0, -17.1), "GG/CC": (-9.3, -23.2), "GT/CA": (-5.9, -12.3),
            "TA/AT": (-7.8, -23.2), "TC/AG": (-5.5, -13.5), "TG/AC": (-9.0, -26.1),
            "TT/AA": (-7.8, -21.9)}

        # Internal mismatch and inosine table (DNA)
        # Allawi & SantaLucia (1997), Biochemistry 36: 10581-10594
        # Allawi & SantaLucia (1998), Biochemistry 37: 9435-9444
        # Allawi & SantaLucia (1998), Biochemistry 37: 2170-2179
        # Allawi & SantaLucia (1998), Nucl Acids Res 26: 2694-2701
        # Peyret et al. (1999), Biochemistry 38: 3468-3477
        # Watkins & SantaLucia (2005), Nucl Acids Res 33: 6258-6267
        DNA_IMM1 = {
            "AG/TT": (1.0, 0.9), "AT/TG": (-2.5, -8.3), "CG/GT": (-4.1, -11.7),
            "CT/GG": (-2.8, -8.0), "GG/CT": (3.3, 10.4), "GG/TT": (5.8, 16.3),
            "GT/CG": (-4.4, -12.3), "GT/TG": (4.1, 9.5), "TG/AT": (-0.1, -1.7),
            "TG/GT": (-1.4, -6.2), "TT/AG": (-1.3, -5.3), "AA/TG": (-0.6, -2.3),
            "AG/TA": (-0.7, -2.3), "CA/GG": (-0.7, -2.3), "CG/GA": (-4.0, -13.2),
            "GA/CG": (-0.6, -1.0), "GG/CA": (0.5, 3.2), "TA/AG": (0.7, 0.7),
            "TG/AA": (3.0, 7.4),
            "AC/TT": (0.7, 0.2), "AT/TC": (-1.2, -6.2), "CC/GT": (-0.8, -4.5),
            "CT/GC": (-1.5, -6.1), "GC/CT": (2.3, 5.4), "GT/CC": (5.2, 13.5),
            "TC/AT": (1.2, 0.7), "TT/AC": (1.0, 0.7),
            "AA/TC": (2.3, 4.6), "AC/TA": (5.3, 14.6), "CA/GC": (1.9, 3.7),
            "CC/GA": (0.6, -0.6), "GA/CC": (5.2, 14.2), "GC/CA": (-0.7, -3.8),
            "TA/AC": (3.4, 8.0), "TC/AA": (7.6, 20.2),
            "AA/TA": (1.2, 1.7), "CA/GA": (-0.9, -4.2), "GA/CA": (-2.9, -9.8),
            "TA/AA": (4.7, 12.9), "AC/TC": (0.0, -4.4), "CC/GC": (-1.5, -7.2),
            "GC/CC": (3.6, 8.9), "TC/AC": (6.1, 16.4), "AG/TG": (-3.1, -9.5),
            "CG/GG": (-4.9, -15.3), "GG/CG": (-6.0, -15.8), "TG/AG": (1.6, 3.6),
            "AT/TT": (-2.7, -10.8), "CT/GT": (-5.0, -15.8), "GT/CT": (-2.2, -8.4),
            "TT/AT": (0.2, -1.5),
            "AI/TC": (-8.9, -25.5), "TI/AC": (-5.9, -17.4), "AC/TI": (-8.8, -25.4),
            "TC/AI": (-4.9, -13.9), "CI/GC": (-5.4, -13.7), "GI/CC": (-6.8, -19.1),
            "CC/GI": (-8.3, -23.8), "GC/CI": (-5.0, -12.6),
            "AI/TA": (-8.3, -25.0), "TI/AA": (-3.4, -11.2), "AA/TI": (-0.7, -2.6),
            "TA/AI": (-1.3, -4.6), "CI/GA": (2.6, 8.9), "GI/CA": (-7.8, -21.1),
            "CA/GI": (-7.0, -20.0), "GA/CI": (-7.6, -20.2),
            "AI/TT": (0.49, -0.7), "TI/AT": (-6.5, -22.0), "AT/TI": (-5.6, -18.7),
            "TT/AI": (-0.8, -4.3), "CI/GT": (-1.0, -2.4), "GI/CT": (-3.5, -10.6),
            "CT/GI": (0.1, -1.0), "GT/CI": (-4.3, -12.1),
            "AI/TG": (-4.9, -15.8), "TI/AG": (-1.9, -8.5), "AG/TI": (0.1, -1.8),
            "TG/AI": (1.0, 1.0), "CI/GG": (7.1, 21.3), "GI/CG": (-1.1, -3.2),
            "CG/GI": (5.8, 16.9), "GG/CI": (-7.6, -22.0),
            "AI/TI": (-3.3, -11.9), "TI/AI": (0.1, -2.3), "CI/GI": (1.3, 3.0),
            "GI/CI": (-0.5, -1.3)}

        # Terminal mismatch table (DNA)
        # SantaLucia & Peyret (2001) Patent Application WO 01/94611
        DNA_TMM1 = {
            "AA/TA": (-3.1, -7.8), "TA/AA": (-2.5, -6.3), "CA/GA": (-4.3, -10.7),
            "GA/CA": (-8.0, -22.5),
            "AC/TC": (-0.1, 0.5), "TC/AC": (-0.7, -1.3), "CC/GC": (-2.1, -5.1),
            "GC/CC": (-3.9, -10.6),
            "AG/TG": (-1.1, -2.1), "TG/AG": (-1.1, -2.7), "CG/GG": (-3.8, -9.5),
            "GG/CG": (-0.7, -19.2),
            "AT/TT": (-2.4, -6.5), "TT/AT": (-3.2, -8.9), "CT/GT": (-6.1, -16.9),
            "GT/CT": (-7.4, -21.2),
            "AA/TC": (-1.6, -4.0), "AC/TA": (-1.8, -3.8), "CA/GC": (-2.6, -5.9),
            "CC/GA": (-2.7, -6.0), "GA/CC": (-5.0, -13.8), "GC/CA": (-3.2, -7.1),
            "TA/AC": (-2.3, -5.9), "TC/AA": (-2.7, -7.0),
            "AC/TT": (-0.9, -1.7), "AT/TC": (-2.3, -6.3), "CC/GT": (-3.2, -8.0),
            "CT/GC": (-3.9, -10.6), "GC/CT": (-4.9, -13.5), "GT/CC": (-3.0, -7.8),
            "TC/AT": (-2.5, -6.3), "TT/AC": (-0.7, -1.2),
            "AA/TG": (-1.9, -4.4), "AG/TA": (-2.5, -5.9), "CA/GG": (-3.9, -9.6),
            "CG/GA": (-6.0, -15.5), "GA/CG": (-4.3, -11.1), "GG/CA": (-4.6, -11.4),
            "TA/AG": (-2.0, -4.7), "TG/AA": (-2.4, -5.8),
            "AG/TT": (-3.2, -8.7), "AT/TG": (-3.5, -9.4), "CG/GT": (-3.8, -9.0),
            "CT/GG": (-6.6, -18.7), "GG/CT": (-5.7, -15.9), "GT/CG": (-5.9, -16.1),
            "TG/AT": (-3.9, -10.5), "TT/AG": (-3.6, -9.8)}

        # Dangling ends table (DNA)
        # Bommarito et al. (2000), Nucl Acids Res 28: 1929-1934
        DNA_DE1 = {
            "AA/.T": (0.2, 2.3), "AC/.G": (-6.3, -17.1), "AG/.C": (-3.7, -10.0),
            "AT/.A": (-2.9, -7.6), "CA/.T": (0.6, 3.3), "CC/.G": (-4.4, -12.6),
            "CG/.C": (-4.0, -11.9), "CT/.A": (-4.1, -13.0), "GA/.T": (-1.1, -1.6),
            "GC/.G": (-5.1, -14.0), "GG/.C": (-3.9, -10.9), "GT/.A": (-4.2, -15.0),
            "TA/.T": (-6.9, -20.0), "TC/.G": (-4.0, -10.9), "TG/.C": (-4.9, -13.8),
            "TT/.A": (-0.2, -0.5),
            ".A/AT": (-0.7, -0.8), ".C/AG": (-2.1, -3.9), ".G/AC": (-5.9, -16.5),
            ".T/AA": (-0.5, -1.1), ".A/CT": (4.4, 14.9), ".C/CG": (-0.2, -0.1),
            ".G/CC": (-2.6, -7.4), ".T/CA": (4.7, 14.2), ".A/GT": (-1.6, -3.6),
            ".C/GG": (-3.9, -11.2), ".G/GC": (-3.2, -10.4), ".T/GA": (-4.1, -13.1),
            ".A/TT": (2.9, 10.4), ".C/TG": (-4.4, -13.1), ".G/TC": (-5.2, -15.0),
            ".T/TA": (-3.8, -12.6)}

        # Dangling ends table (RNA)
        # Turner & Mathews (2010), Nucl Acids Res 38: D280-D282
        RNA_DE1 = {
            ".T/AA": (-4.9, -13.2), ".T/CA": (-0.9, -1.3), ".T/GA": (-5.5, -15.1),
            ".T/TA": (-2.3, -5.5),
            ".G/AC": (-9.0, -23.5), ".G/CC": (-4.1, -10.6), ".G/GC": (-8.6, -22.2),
            ".G/TC": (-7.5, -20.31),
            ".C/AG": (-7.4, -20.3), ".C/CG": (-2.8, -7.7), ".C/GG": (-6.4, -16.4),
            ".C/TG": (-3.6, -9.7),
            ".T/AG": (-4.9, -13.2), ".T/CG": (-0.9, -1.3), ".T/GG": (-5.5, -15.1),
            ".T/TG": (-2.3, -5.5),
            ".A/AT": (-5.7, -16.1), ".A/CT": (-0.7, -1.9), ".A/GT": (-5.8, -16.4),
            ".A/TT": (-2.2, -6.8),
            ".G/AT": (-5.7, -16.1), ".G/CT": (-0.7, -1.9), ".G/GT": (-5.8, -16.4),
            ".G/TT": (-2.2, -6.8),
            "AT/.A": (-0.5, -0.6), "CT/.A": (6.9, 22.6), "GT/.A": (0.6, 2.6),
            "TT/.A": (0.6, 2.6),
            "AG/.C": (-1.6, -4.5), "CG/.C": (0.7, 3.2), "GG/.C": (-4.6, -14.8),
            "TG/.C": (-0.4, -1.3),
            "AC/.G": (-2.4, -6.1), "CC/.G": (3.3, 11.6), "GC/.G": (0.8, 3.2),
            "TC/.G": (-1.4, -4.2),
            "AT/.G": (-0.5, -0.6), "CT/.G": (6.9, 22.6), "GT/.G": (0.6, 2.6),
            "TT/.G": (0.6, 2.6),
            "AA/.T": (1.6, 6.1), "CA/.T": (2.2, 8.1), "GA/.T": (0.7, 3.5),
            "TA/.T": (3.1, 10.6),
            "AG/.T": (1.6, 6.1), "CG/.T": (2.2, 8.1), "GG/.T": (0.7, 3.5),
            "TG/.T": (3.1, 10.6)}
        # Set defaults
        if not nn_table:
            nn_table = DNA_NN3
        if not tmm_table:
            tmm_table = DNA_TMM1
        if not imm_table:
            imm_table = DNA_IMM1
        if not de_table:
            de_table = DNA_DE1

        c_seq = Essentials.Antisense(seq)

        tmp_seq = seq
        tmp_cseq = c_seq
        delta_h = 0
        delta_s = 0
        d_h = 0  # Names for indexes
        d_s = 1  # 0 and 1

        # Now for terminal mismatches
        left_tmm = tmp_cseq[:2][::-1] + "/" + tmp_seq[:2][::-1]
        if left_tmm in tmm_table:
            delta_h += tmm_table[left_tmm][d_h]
            delta_s += tmm_table[left_tmm][d_s]
            tmp_seq = tmp_seq[1:]
            tmp_cseq = tmp_cseq[1:]
        right_tmm = tmp_seq[-2:] + "/" + tmp_cseq[-2:]
        if right_tmm in tmm_table:
            delta_h += tmm_table[right_tmm][d_h]
            delta_s += tmm_table[right_tmm][d_s]
            tmp_seq = tmp_seq[:-1]
            tmp_cseq = tmp_cseq[:-1]

        # Now everything 'unusual' at the ends is handled and removed and we can
        # look at the initiation.
        # One or several of the following initiation types may apply:

        # Type: General initiation value
        delta_h += nn_table["init"][d_h]
        delta_s += nn_table["init"][d_s]

        # Type: Duplex with no (allA/T) or at least one (oneG/C) GC pair
        if GC == 0:
            delta_h += nn_table["init_allA/T"][d_h]
            delta_s += nn_table["init_allA/T"][d_s]
        else:
            delta_h += nn_table["init_oneG/C"][d_h]
            delta_s += nn_table["init_oneG/C"][d_s]

        # Type: Penalty if 5' end is T
        if seq.startswith("T"):
            delta_h += nn_table["init_5T/A"][d_h]
            delta_s += nn_table["init_5T/A"][d_s]
        if seq.endswith("A"):
            delta_h += nn_table["init_5T/A"][d_h]
            delta_s += nn_table["init_5T/A"][d_s]

        # Type: Different values for G/C or A/T terminal basepairs
        ends = seq[0] + seq[-1]
        AT = ends.count("A") + ends.count("T")
        GC = ends.count("G") + ends.count("C")
        delta_h += nn_table["init_A/T"][d_h] * AT
        delta_s += nn_table["init_A/T"][d_s] * AT
        delta_h += nn_table["init_G/C"][d_h] * GC
        delta_s += nn_table["init_G/C"][d_s] * GC

        # Finally, the 'zipping'
        for basenumber in range(len(tmp_seq) - 1):
            neighbors = (
                tmp_seq[basenumber : basenumber + 2]
                + "/"
                + tmp_cseq[basenumber : basenumber + 2]
            )
            if neighbors in imm_table:
                delta_h += imm_table[neighbors][d_h]
                delta_s += imm_table[neighbors][d_s]
            elif neighbors[::-1] in imm_table:
                delta_h += imm_table[neighbors[::-1]][d_h]
                delta_s += imm_table[neighbors[::-1]][d_s]
            elif neighbors in nn_table:
                delta_h += nn_table[neighbors][d_h]
                delta_s += nn_table[neighbors][d_s]
            elif neighbors[::-1] in nn_table:
                delta_h += nn_table[neighbors[::-1]][d_h]
                delta_s += nn_table[neighbors[::-1]][d_s]
            else:
                # We haven't found the key...
                _key_error(neighbors, strict)

        k = (dnac1 - (dnac2 / 2.0)) * 1e-9
        R = 1.987  # universal gas constant in Cal/degrees C*Mol
        if salt_correction:
            def get_salt_correction(Na = Na,
                                    K = K,
                                    Tris = Tris,
                                    Mg = Mg,
                                    dNTPs = dNTPs):
                Na=50
                K=0
                Tris=0
                Mg=0
                dNTPs = 0

                Mon = Na + K + Tris / 2.0 
                mon = Mon * 1e-3
                corr = 16.6 * np.log10(mon)

                return corr

            corr = get_salt_correction()
            delta_s += corr

        melting_temp = (1000 * delta_h) / (delta_s + (R * (np.log(k)))) - 273.15

        return melting_temp

    @jit(nopython=True)
    def increment(num, maximum):
        num = num.copy()
        num[-1]+=1
        if num[-1]>maximum[-1]:
            for digit in range(len(num)):
                if num[-digit]>maximum[-digit]:
                    num[-digit] = 0
                    num[-digit-1] += 1
        return num.copy()

    def get_all_possible_versions(seq):       
        mapping = Essentials.IUPAC.copy()

        if "N" in seq:
            for key in mapping.keys():
                mapping[key].append("N")
                
        maximum = np.zeros(len(seq))
        for i in range(len(seq)):
            maximum[i] = len(mapping[seq[i]])-1
        
        encoded = [np.zeros(len(seq))]
        while (encoded[-1]!=maximum).any():
            encoded.append(Essentials.increment(encoded[-1], maximum))
        
        res = []
        for enc in encoded:
            primer = ""
            for offset, i in enumerate(enc):
                primer+=mapping[seq[offset]][int(i)]
            res.append(primer)
        return res
