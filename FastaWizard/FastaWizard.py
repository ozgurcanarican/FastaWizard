#__author__: Özgür Can ARICAN
#__contact__: ozgurcanarican95@gmail.com

"""
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import time
import sys
import os

print("""FastaWizard v_1.0  Copyright © 2020  Özgür Can ARICAN
This program comes with ABSOLUTELY NO WARRANTY; for details type `w'.
This is free software, and you are welcome to redistribute it under 
certain conditions; type `c' for details.
Type 'a' for detail about program and contact information about author""")

os.system("color")

print("""\033[34m
***************************************************************
_______                    _       _                         _
| _____|         _        | |     | |                       | |
| |___ ___ _ ___| |_ ___ _| | / \ | |_ ____  ___ _ _ __  ___| |
| ____/ _ ' / __| __/ _ ' | |/   \| | |__  |/ _ ' | '__|/ _   |
| |  | (_|  \__ \ |_ (_|  | / / \ \ | | / /  (_|  | |  | (_|  |
|_|   \___._|___/\__\___._|__/   \__|_|/___|\___._|_|   \___._|

***************************************************************
""")

initiater = input("Press 'Enter' to start program...\033[0m")
if(initiater == "w"):
    print("""\nFastaWizard v_1.0  Copyright © 2020  Özgür Can ARICAN
    
    Disclaimer of Warranty.

THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY
APPLICABLE LAW.  EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT
HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM "AS IS" WITHOUT WARRANTY
OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE.  THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM
IS WITH YOU.  SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF
ALL NECESSARY SERVICING, REPAIR OR CORRECTION.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
    """)
    input("Press 'Enter' to start program...")
elif(initiater == "c"):
    print("""\nFastaWizard v_1.0  Copyright © 2020  Özgür Can ARICAN
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
    """)
    input("Press 'Enter' to start program...")
elif(initiater == "a"):
    print("""
    -------------------------------------------------
    |   Program Name     |       FastaWizard        |
    |      Version       |           1.0            |
    |   Completion Date  |       2nd Jan 2020       |
    |    Developed at    |        Python 3.8        |
    |   Author's Name    |      Özgür Can ARICAN    |
    |    Gmail Adress    |ozgurcanarican95@gmail.com|
    |      GitHub        |github.com/ozgurcanarican |
    -------------------------------------------------
    """)
    input("Press 'Enter' to start program...")

class FastaWizard:

    def __init__(self):
        self.sequence = ""
        self.sequenceTwo = ""
        self.complementary = ""
        self.codon = ""
        self.anticodon = ""
        self.sequence_triplet = []
        self.triple_aa_list = []
        self.single_aa_list = []
        self.detect_triple_list = []
        self.detect_single_list = []
        self.numerator = 0
        self.denominator = 0
        self.similarity = 0

    @staticmethod
    def template_order():
        return sequence.upper()

    def template_to_complementary(self):
        self.sequence_comp = sequence
        for i in self.sequence_comp:
            if (i == "a" or i == "A"):
                i = "T"
            elif (i == "t" or i == "T"):
                i = "A"
            elif (i == "g" or i == "G"):
                i = "C"
            elif (i == "c" or i == "C"):
                i = "G"
            else:
                pass
            self.complementary = self.complementary + i
        return self.complementary

    def template_to_codon(self):
        self.sequence_cod = sequence
        for j in self.sequence_cod:
            if (j == "a" or j == "A"):
                j = "U"
            elif (j == "t" or j == "T"):
                j = "A"
            elif (j == "g" or j == "G"):
                j = "C"
            elif (j == "c" or j == "C"):
                j = "G"
            else:
                pass
            self.codon = self.codon + j
        return self.codon

    def template_to_anticodon(self):
        self.sequence_anti = sequence
        for k in self.sequence_anti:
            if (k == "a" or k == "A"):
                k = "A"
            elif (k == "t" or k == "T"):
                k = "U"
            elif (k == "g" or k == "G"):
                k = "G"
            elif (k == "c" or k == "C"):
                k = "C"
            else:
                pass
            self.anticodon = self.anticodon + k
        return self.anticodon

    @property
    def template_to_aminoacid(self):
        self.template_to_anticodon()
        self.sequence_triplet = [self.anticodon[x:x + 3] for x in range(0, len(self.anticodon), 3)]

        self.aa_dict_three_letters = {"UUU": "Phe", "UUC": "Phe", "UUA": "Leu", "UUG": "Leu", "UCU": "Ser", "UCC": "Ser",
                                      "UCA": "Ser", "UCG": "Ser","UAU": "Tyr", "UAC": "Tyr", "UAA": "Stop", "UAG": "Stop",
                                      "UGU": "Cys", "UGC": "Cys", "UGA": "Stop", "UGG": "Trp", "CUU": "Leu", "CUC": "Leu",
                                      "CUA": "Leu", "CUG": "Leu", "CCU": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
                                      "CAU": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln", "CGU": "Arg", "CGC": "Arg",
                                      "CGA": "Arg", "CGG": "Arg", "AUU": "Ile", "AUC": "Ile", "AUA": "Ile", "AUG": "Met",
                                      "ACU": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr", "AAU": "Asn", "AAC": "Asn",
                                      "AAA": "Lys", "AAG": "Lys", "AGU": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg",
                                      "GUU": "Val", "GUC": "Val", "GUA": "Val", "GUG": "Val", "GCU": "Ala", "GCC": "Ala",
                                      "GCA": "Ala", "GCG": "Ala", "GAU": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
                                      "GGU": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly"}

        self.aa_dict_single_letter = {"UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L", "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
                                      "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*", "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
                                      "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L", "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
                                      "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
                                      "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M", "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
                                      "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
                                      "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V", "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
                                      "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E", "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"}

        for firstVar in range(len(self.sequence_triplet)):
            if self.sequence_triplet[firstVar] in self.aa_dict_three_letters:
                self.triple_aa_list.append(self.aa_dict_three_letters[self.sequence_triplet[firstVar]])
            else:
                self.triple_aa_list.append('NoNe')

        for secondVar in range(len(self.sequence_triplet)):
            if self.sequence_triplet[secondVar] in self.aa_dict_single_letter:
                self.single_aa_list.append(self.aa_dict_single_letter[self.sequence_triplet[secondVar]])
            else:
                self.single_aa_list.append('.')

        self.start_stop_codon()
        for thirdVar in range(len(self.sequence_triplet)):
            if self.sequence_triplet[thirdVar] in self.aa_dict_three_letters:
                self.detect_triple_list.append(self.aa_dict_three_letters[self.sequence_triplet[thirdVar]])
            else:
                self.detect_triple_list.append('NoNe')

        for forthVar in range(len(self.sequence_triplet)):
            if self.sequence_triplet[forthVar] in self.aa_dict_single_letter:
                self.detect_single_list.append(self.aa_dict_single_letter[self.sequence_triplet[forthVar]])
            else:
                self.detect_single_list.append('.')

        return "-".join(self.triple_aa_list), "".join(self.single_aa_list), "-".join(self.detect_triple_list), "".join(self.detect_single_list)

    def start_stop_codon(self):
        for fifthVar in self.sequence_triplet:
            self.sequence_triplet.remove(fifthVar)
            if fifthVar == "AUG":
                for sixthVar in self.sequence_triplet:
                    if sixthVar == "UAA" or sixthVar == "UAG" or sixthVar == "UGA":
                        break

    @property
    def adenine_ratio(self):
        self.sequence_adenine = sequence
        for a in self.sequence_adenine:
            if (a == "a" or a == "A"):
                self.numerator += 1
            self.denominator += 1
        self.adenineRatio = (self.numerator / self.denominator) * 100
        self.numerator = 0
        self.denominator = 0
        return self.adenineRatio

    @property
    def thymine_ratio(self):
        self.sequence_thymine = sequence
        for t in self.sequence_thymine:
            if (t == "t" or t == "T"):
                self.numerator += 1
            self.denominator += 1
        self.thymineRatio = (self.numerator / self.denominator) * 100
        self.numerator = 0
        self.denominator = 0
        return self.thymineRatio

    @property
    def guanine_ratio(self):
        self.sequence_guanine = sequence
        for g in self.sequence_guanine:
            if (g == "g" or g == "G"):
                self.numerator += 1
            self.denominator += 1
        self.guanineRatio = (self.numerator / self.denominator) * 100
        self.numerator = 0
        self.denominator = 0
        return self.guanineRatio

    @property
    def cytosine_ratio(self):
        self.sequence_cytosine = sequence
        for c in self.sequence_cytosine:
            if (c == "c" or c == "C"):
                self.numerator += 1
            self.denominator += 1
        self.cytosineRatio = (self.numerator / self.denominator) * 100
        self.numerator = 0
        self.denominator = 0
        return self.cytosineRatio

    @property
    def uracile_ratio(self):
        self.sequence_uracile = sequence
        for u in self.sequence_uracile:
            if (u == "u" or u == "U"):
                self.numerator += 1
            self.denominator += 1
        self.uracileRatio = (self.numerator / self.denominator) * 100
        self.numerator = 0
        self.denominator = 0
        return self.uracileRatio

    @property
    def gc_content(self):
        self.sequence_gc = sequence
        for gc in self.sequence_gc:
            if (gc == "g" or gc == "G" or gc == "c" or gc == "C"):
                self.numerator += 1
            self.denominator += 1
        self.gcContent = (self.numerator / self.denominator) * 100
        self.numerator = 0
        self.denominator = 0
        return self.gcContent

    @property
    def at_content(self):
        self.sequence_at = sequence
        for at in self.sequence_at:
            if (at == "a" or at == "A" or at == "t" or at == "T"):
                self.numerator += 1
            self.denominator += 1
        self.atContent = (self.numerator / self.denominator) * 100
        self.numerator = 0
        self.denominator = 0
        return self.atContent

    def sequence_similarity(self):
        self.zipped = zip(sequence, sequenceTwo)
        for n1, n2 in self.zipped:
            if (n1 == n2):
                self.numerator += 1
            self.denominator += 1
        self.similarity = (self.numerator / self.denominator) * 100
        self.numerator = 0
        self.denominator = 0
        return self.similarity


def sequence_options():
    global sequence
    global header
    global sequenceTwo

    if (process == "1" or process == "2"):
        print("\n\n\033[36mATTENTION: if you enter a fasta file, your file must be at the same directory with the program!!!\033[0m")
        while True:
            isFasta = input("\nDo you enter your sequence as a fasta file ? (y/n): ")
            if (isFasta == "n" or isFasta == "no" or isFasta == "N" or isFasta == "No" or isFasta == "NO" or isFasta == "y" or isFasta == "yes" or isFasta == "Y" or isFasta == "Yes" or isFasta == "YES"):
                time.sleep(1)
                break
            else:
                print("\033[31mERROR: Invalid answer. Please enter your answer as y, yes, Y, Yes, YES, n, no, N, No or NO...\033[0m")
                time.sleep(2)
                continue
        while True:
            cutSequence = input("\nWould you like to specify the location of your sequence for the process ? (y/n): ")
            if (cutSequence == "n" or cutSequence == "no" or cutSequence == "N" or cutSequence == "No" or cutSequence == "NO" or cutSequence == "y" or cutSequence == "yes" or cutSequence == "Y" or cutSequence == "Yes" or cutSequence == "YES"):
                time.sleep(1)
                break
            else:
                print("\033[31mERROR: Invalid answer. Please enter your answer as y, yes, Y, Yes, YES, n, no, N, No or NO...\033[0m")
                time.sleep(2)
                continue

        if ((isFasta == "n" or isFasta == "no" or isFasta == "N" or isFasta == "No" or isFasta == "NO") and (cutSequence == "n" or cutSequence == "no" or cutSequence == "N" or cutSequence == "No" or cutSequence == "NO")):
            validationOne = True
            while validationOne:
                sequence = input("\nPlease enter manually or copy/paste your sequence: ")
                for valueOne in sequence:
                    if (valueOne == "a" or valueOne == "A" or valueOne == "t" or valueOne == "T" or valueOne == "g" or valueOne == "G" or valueOne == "c" or valueOne == "C" or valueOne == "u" or valueOne == "U"):
                        pass
                    else:
                        print("\033[31mERROR: Invalid sequence entry. Please check your sequence if there are any spaces, integers or cursors\033[0m")
                        time.sleep(2)
                        break
                else:
                    time.sleep(1)
                    validationOne = False

        elif ((isFasta == "y" or isFasta == "yes" or isFasta == "Y" or isFasta == "Yes" or isFasta == "YES") and (cutSequence == "n" or cutSequence == "no" or cutSequence == "N" or cutSequence == "No" or cutSequence == "NO")):
            while True:
                sequence_entry = input("\nPlease enter your fasta file name, with the file extension (ex. sequence.fasta): ")
                if (sequence_entry.endswith(".fasta") == True):
                    try:
                        sequence_file = open(sequence_entry, 'r')
                        time.sleep(1)
                        break
                    except FileNotFoundError:
                        print("\033[31mERROR: Please check if the fasta file name correct or the file is at the same directory with the program...\033[0m")
                        time.sleep(2)
                        continue
                else:
                    print("\033[31mERROR: Please check if the fasta file extension is correct. It should ends with '.fasta'\033[0m")
                    time.sleep(2)
                    continue
            sequence_parse = sequence_file.readline()
            while sequence_parse:
                sequence_parse = sequence_parse.rstrip("\n")
                if ">" in sequence_parse:
                    header = sequence_parse
                else:
                    sequence = sequence + sequence_parse
                sequence_parse = sequence_file.readline()

        elif ((isFasta == "n" or isFasta == "no" or isFasta == "N" or isFasta == "No" or isFasta == "NO") and (cutSequence == "y" or cutSequence == "yes" or cutSequence == "Y" or cutSequence == "Yes" or cutSequence == "YES")):
            validationTwo = True
            while validationTwo:
                sequence = input("\nPlease enter manually or copy/paste your sequence: ")
                for valueTwo in sequence:
                    if (valueTwo == "a" or valueTwo == "A" or valueTwo == "t" or valueTwo == "T" or valueTwo == "g" or valueTwo == "G" or valueTwo == "c" or valueTwo == "C" or valueTwo == "u" or valueTwo == "U"):
                        pass
                    else:
                        print("ERROR: Invalid sequence entry. Please check your sequence if there are any spaces, integers or cursors")
                        time.sleep(2)
                        break
                else:
                    time.sleep(1)
                    validationTwo = False
            while True:
                try:
                    pos1 = int(input("\ncut dna from position: "))
                    time.sleep(1)
                    break
                except ValueError:
                    print("\033[31mERROR: Please check your input. It should be an integer\033[0m")
                    time.sleep(2)
                    continue
            while True:
                try:
                    pos2 = int(input("to position: "))
                    time.sleep(1)
                    break
                except ValueError:
                    print("\33[31mERROR: Please check your input. It should be an integer\033[0m")
                    time.sleep(2)
                    continue
            modSequence = ""
            counter = 0
            for nucleotide in sequence:
                counter += 1
                if (pos1 <= counter <= pos2):
                    modSequence = modSequence + nucleotide
            sequence = modSequence

        elif ((isFasta == "y" or isFasta == "yes" or isFasta == "Y" or isFasta == "Yes" or isFasta == "YES") and (cutSequence == "y" or cutSequence == "yes" or cutSequence == "Y" or cutSequence == "Yes" or cutSequence == "YES")):
            sequence = ""
            header = ""
            while True:
                sequence_entry = input("\nPlease enter your fasta file name, with the file extension (ex. sequence.fasta): ")
                if (sequence_entry.endswith(".fasta") == True):
                    try:
                        sequence_file = open(sequence_entry, 'r')
                        break
                    except FileNotFoundError:
                        print("\033[31mERROR: Please check if the fasta file name correct or the file is at the same directory with the program...\033[0m")
                        time.sleep(4)
                        continue
                else:
                    print("\033[31mERROR: Please check if the fasta file extension is correct. It should ends with '.fasta'\033[0m")
                    time.sleep(4)
                    continue
            sequence_parse = sequence_file.readline()
            while True:
                try:
                    pos1 = int(input("\ncut dna from position: "))
                    time.sleep(1)
                    break
                except ValueError:
                    print("\033[31mERROR: Please check your input. It should be an integer\033[0m")
                    time.sleep(2)
                    continue
            while True:
                try:
                    pos2 = int(input("to position: "))
                    time.sleep(1)
                    break
                except ValueError:
                    print("\33[31mERROR: Please check your input. It should be an integer\033[0m")
                    time.sleep(2)
                    continue
            while sequence_parse:
                sequence_parse = sequence_parse.rstrip("\n")
                if ">" in sequence_parse:
                    header = sequence_parse
                else:
                    sequence = sequence + sequence_parse
                sequence_parse = sequence_file.readline()
                modSequence = ""
                counter = 0
                for nucleotide in sequence:
                    counter += 1
                    if (pos1 <= counter <= pos2):
                        modSequence = modSequence + nucleotide
                sequence = modSequence


    elif (process == "3"):
        print("\n\n\033[36mNote: if you enter fasta files, your files must be at the same directory with the program!!!\033[0m")
        while True:
            isFasta = input("\nDo you enter your sequences as fasta files ? (y/n): ")
            if (isFasta == "n" or isFasta == "no" or isFasta == "N" or isFasta == "No" or isFasta == "NO" or isFasta == "y" or isFasta == "yes" or isFasta == "Y" or isFasta == "Yes" or isFasta == "YES"):
                time.sleep(1)
                break
            else:
                print("\033[31mERROR: Invalid answer. Please enter your answer as y, yes, Y, Yes, YES, n, no, N, No or NO...\033[0m")
                time.sleep(2)
                continue

        if (isFasta == "n" or isFasta == "no" or isFasta == "N" or isFasta == "No" or isFasta == "NO"):
            validationThree = True
            while validationThree:
                sequence = input("\nPlease enter manually or copy/paste your first sequence: ")
                for valueThree in sequence:
                    if (valueThree == "a" or valueThree == "A" or valueThree == "t" or valueThree == "T" or valueThree == "g" or valueThree == "G" or valueThree == "c" or valueThree == "C" or valueThree == "u" or valueThree == "U"):
                        pass
                    else:
                        print("\033[31mERROR: Invalid sequence entry. Please check your sequence if there are any spaces, integers or cursors\033[0m")
                        time.sleep(2)
                        break
                else:
                    time.sleep(1)
                    validationThree = False
            validationFour = True
            while validationFour:
                sequenceTwo = input("Please enter the second sequence: ")
                for valueFour in sequence:
                    if (valueFour == "a" or valueFour == "A" or valueFour == "t" or valueFour == "T" or valueFour == "g" or valueFour == "G" or valueFour == "c" or valueFour == "C" or valueFour == "u" or valueFour == "U"):
                        pass
                    else:
                        print("\033[31mERROR: Invalid sequence entry. Please check your sequence if there are any spaces, integers or cursors\033[0m")
                        time.sleep(2)
                        break
                else:
                    time.sleep(1)
                    validationFour = False


        elif (isFasta == "y" or isFasta == "yes" or isFasta == "Y" or isFasta == "Yes" or isFasta == "YES"):
            while True:
                sequence_entry = input("\nPlease enter your fasta file name, with the file extension (ex. sequence.fasta): ")
                if (sequence_entry.endswith(".fasta") == True):
                    try:
                        sequence_file = open(sequence_entry, 'r')
                        time.sleep(1)
                        break
                    except FileNotFoundError:
                        print("\033[31mERROR: Please check if the fasta file name correct or the file is at the same directory with the program...\033[0m")
                        time.sleep(2)
                        continue
            sequence_parse = sequence_file.readline()
            while sequence_parse:
                sequence_parse = sequence_parse.rstrip("\n")
                if ">" in sequence_parse:
                    header = sequence_parse
                else:
                    sequence = sequence + sequence_parse
                sequence_parse = sequence_file.readline()
            header = ""
            while True:
                sequenceTwo_entry = input("\nPlease enter the other fasta file name with the extension: ")
                if (sequenceTwo_entry.endswith(".fasta") == True):
                    try:
                        sequenceTwo_file = open(sequenceTwo_entry, 'r')
                        time.sleep(1)
                        break
                    except FileNotFoundError:
                        print("\033[31mERROR: Please check if the fasta file name correct or the file is at the same directory with the program...\033[0m")
                        time.sleep(2)
                        continue
            sequenceTwo_parse = sequenceTwo_file.readline()
            while sequenceTwo_parse:
                sequenceTwo_parse = sequenceTwo_parse.rstrip("\n")
                if ">" in sequenceTwo_parse:
                    header = sequenceTwo_parse
                else:
                    sequenceTwo = sequenceTwo + sequenceTwo_parse
                sequenceTwo_parse = sequenceTwo_file.readline()


def save_options():
    while True:
        saveFile = input("\n\nWould you like to save the results as text file ? (y/n): ")
        if (saveFile == "n" or saveFile == "no" or saveFile == "N" or saveFile == "No" or saveFile == "NO" or saveFile == "y" or saveFile == "yes" or saveFile == "Y" or saveFile == "Yes" or saveFile == "YES"):
            time.sleep(1)
            break
        else:
            print("\033[31mERROR: Invalid answer. Please enter your answer as y, yes, Y, Yes, YES, n, no, N, No or NO...\033[0m")
            time.sleep(2)
            continue
    if (saveFile == "n" or saveFile == "no" or saveFile == "N" or saveFile == "No" or saveFile == "NO"):
        time.sleep(3)
        pass
    elif (saveFile == "y" or saveFile == "yes" or saveFile == "Y" or saveFile == "Yes" or saveFile == "YES"):
        time.sleep(3)
        if (subprocessOne == "1"):
            with open("complementary.txt", "a") as text_file:
                text_file.write("The complementary of your given sequence is:\n" + fastaWizard.template_to_complementary())
        elif (subprocessOne == "2"):
            with open("mRNA_translation.txt", "a") as text_file:
                text_file.write("The codon (mRNA sequence) of your given sequence is:\n" + fastaWizard.template_to_codon())
        elif (subprocessOne == "3"):
            with open("tRNA_translation.txt", "a") as text_file:
                text_file.write("The anticodon (tRNA sequence) of your given sequence is:\n" + fastaWizard.template_to_anticodon())
        elif (subprocessOne == "4" and aaType == "1" and (codonDetection == "n" or codonDetection == "N" or codonDetection == "no" or codonDetection == "No" or codonDetection == "NO")):
            with open("single-letter_aa_translation.txt", "a") as text_file:
                text_file.write("The single-letter aminoacid sequence of your given sequence is:\n" + fastaWizard.template_to_aminoacid[1])
        elif (subprocessOne == "4" and aaType == "1" and (codonDetection == "y" or codonDetection == "Y" or codonDetection == "yes" or codonDetection == "Yes" or codonDetection == "YES")):
            with open("single-letter_aa_translation.txt", "a") as text_file:
                text_file.write("The single-letter aminoacid sequence of your given sequence is:\n" + fastaWizard.template_to_aminoacid[3])
        elif (subprocessOne == "4" and aaType == "2" and (codonDetection == "n" or codonDetection == "N" or codonDetection == "no" or codonDetection == "No" or codonDetection == "NO")):
            with open("three-letters_aa_translation.txt", "a") as text_file:
                text_file.write("The three-letter aminoacid sequence of your given sequence is:\n" + fastaWizard.template_to_aminoacid[0])
        elif (subprocessOne == "4" and aaType == "2" and (codonDetection == "y" or codonDetection == "Y" or codonDetection == "yes" or codonDetection == "Yes" or codonDetection == "YES")):
            with open("three-letters_aa_translation.txt", "a") as text_file:
                text_file.write("The three-letter aminoacid sequence of your given sequence is:\n" + fastaWizard.template_to_aminoacid[2])
        elif (subprocessOne == "4" and aaType == "3" and (codonDetection == "n" or codonDetection == "N" or codonDetection == "no" or codonDetection == "No" or codonDetection == "NO")):
            with open("aa_translations.txt", "a") as text_file:
                text_file.write("The aminoacid sequences of your given sequence are: :\n" + fastaWizard.template_to_aminoacid[1] + "\n" + fastaWizard.template_to_aminoacid[0])
        elif (subprocessOne == "4" and aaType == "3" and (codonDetection == "y" or codonDetection == "Y" or codonDetection == "yes" or codonDetection == "Yes" or codonDetection == "YES")):
            with open("aa_translations.txt", "a") as text_file:
                text_file.write("The aminoacid sequences of your given sequence are: :\n" + fastaWizard.template_to_aminoacid[3] + "\n" + fastaWizard.template_to_aminoacid[2])
        elif (subprocessOne == "5"):
            with open("full_translations.txt", "a") as text_file:
                 text_file.write("The complementary of your given sequence is:\n" + res_all_complementary + "\n\n"
                                 + "The codon (mRNA sequence) of your given sequence is:\n" + res_all_anticodon + "\n\n"
                                 + "The anticodon (tRNA sequence) of your given sequence is:\n" + res_all_anticodon + "\n\n"
                                 + "The single-letter aminoacid sequence of your given sequence is:\n" + res_all_aa_single + "\n\n"
                                 + "The three-letters aminoacid sequence of your given sequence is:\n" + res_all_aa_triple)
        elif (subprocessTwo == "1"):
            with open("g-c_content.txt", "a") as text_file:
                text_file.write("G-C concentration is " + str(res_gc_content) + "% of your given sequence")
        elif (subprocessTwo == "2"):
            with open("a-t_content.txt", "a") as text_file:
                text_file.write("A-T concentration is " + str(res_at_content) + "% of your given sequence")
        elif (subprocessTwo == "3"):
            with open("atgc_ratios.txt", "a") as text_file:
                text_file.write("Adenine ratio = " + "%" + str(res_dna_aratio) + "\n" +
                                "Thymine ratio = " + "%" + str(res_dna_tratio) + "\n" +
                                "Guanine ratio = " + "%" + str(res_dna_gratio) + "\n" +
                                "Cytosine ratio = " + "%" + str(res_dna_cratio))
        elif (subprocessTwo == "4"):
            with open("augc_ratios.txt", "a") as text_file:
                text_file.write("Adenine ratio = " + "%" + str(res_rna_aratio) + "\n" +
                                "Uracile ratio = " + "%" + str(res_rna_uratio) + "\n" +
                                "Guanine ratio = " + "%" + str(res_rna_gratio) + "\n" +
                                "Cytosine ratio = " + "%" + str(res_rna_cratio))
        elif (subprocessTwo == "5"):
            with open("full_content_analysis.txt", "a") as text_file:
                text_file.write("Adenine ratio = " + "%" + str(res_all_aratio) + "\n" +
                                "Thymine ratio = " + "%" + str(res_all_tratio) + "\n" +
                                "Guanine ratio = " + "%" + str(res_all_gratio) + "\n" +
                                "Cytosine ratio = " + "%" + str(res_all_cratio) + "\n\n" +
                                "G-C concentration is " + str(res_all_gc) + "% of your given sequence" + "\n" +
                                "A-T concentration is " + str(res_all_at) + "% of your given sequence")
        elif (process == "3"):
            with open("nucleotide_similarity.txt", "a") as text_file:
                text_file.write("ooooooooooooooooooooooooooooooooooooo" + "\n" +
                                "|||||||||||||||||||||||||||||||||||||" + "\n" +
                                "ooooooooooooooooooooooooooooooooooooo" + "\n\n" +
                                "The nucleotide similarity between your sequences are " + str(similarity) + " %")

        print("\n\nYour text file has generated")


print("\n\033[32mWelcome to FastaWizard. This program is a bioinformatics tool\nthat allows you to analyze on a DNA sequence.\033[0m")
while True:

    print("""
\033[32m**************************************************************************

1) Sequence translation (from template to complementary, codon, anticodon or amino acid sequence)

2) Content analysis (G-C content, A-T content or A, T, G, C, U ratios)

3) Sequence similarity (analyze the similarity between to dna or rna sequences)

Please press 'q' to quit program

**************************************************************************\033[0m\n\n
    """)

    fastaWizard = FastaWizard()
    sequence = ""
    sequenceTwo = ""
    header = ""
    process = input("Please select the process you would like to do: ")

    if (process == "1"):
        sequence_options()
        print("""
\033[32m**************************************
1) Template to complementary

2) Template to codon

3) Template to anticodon

4) Template to amino acid sequence
            
5) All translations above (all options are default)

Please press 'q' to return main menu
**************************************\033[0m\n\n
        """)
        time.sleep(1)
        subprocessOne = input("Please select the subprocess you would like to do: ")

        if (subprocessOne == "1"):
            res_complementary = fastaWizard.template_to_complementary()
            print("\nThe complementary of your given sequence is: ")
            print(res_complementary)
            time.sleep(5)
            save_options()

        elif (subprocessOne == "2"):
            res_codon = fastaWizard.template_to_codon()
            print("\nThe codon (mRNA sequence) of your given sequence is: ")
            print(res_codon)
            time.sleep(5)
            save_options()

        elif (subprocessOne == "3"):
            res_anticodon = fastaWizard.template_to_anticodon()
            print("\nThe anticodon (tRNA sequence) of your given sequence is: ")
            print(res_anticodon)
            time.sleep(5)
            save_options()

        elif (subprocessOne == "4"):
            print("\nWould you prefer the amino acid sequence as single-letter or three-letters ?")
            print("1) Single-letter\n2) Three-letters\n3) Both single and three letters ")
            while True:
                aaType = input("Please select one option: ")
                if (aaType == "1" or aaType == "2" or aaType == "3"):
                    time.sleep(1)
                    break
                else:
                    print("\033[31mERROR: Please choose one option... Do not enter unordinary strings or integers to input box...\033[0m")
                    time.sleep(2)
                    continue
            while True:
                codonDetection = input("\nWould you like the detection of start and stop codon in your sequence ? (y/n): ")
                if (codonDetection == "n" or codonDetection == "no" or codonDetection == "N" or codonDetection == "No" or codonDetection == "NO" or codonDetection == "y" or codonDetection == "yes" or codonDetection == "Y" or codonDetection == "Yes" or codonDetection == "YES"):
                    time.sleep(1)
                    break
                else:
                    print("\033[31mERROR: Invalid answer. Please enter your answer as y, yes, Y, Yes, YES, n, no, N, No or NO...\033[0m")
                    time.sleep(2)
                    continue
            if (aaType == "1" and (codonDetection == "n" or codonDetection == "N" or codonDetection == "no" or codonDetection == "No" or codonDetection == "NO")):
                res_aa_one = fastaWizard.template_to_aminoacid[1]
                print("\nThe single-letter amino acid sequence of your given sequence is: ")
                print(res_aa_one)
                save_options()

            elif (aaType == "1" and (codonDetection == "y" or codonDetection == "Y" or codonDetection == "yes" or codonDetection == "Yes" or codonDetection == "YES")):
                res_aa_three = fastaWizard.template_to_aminoacid[3]
                print("\nThe single-letter amino acid sequence of your given sequence is: ")
                print(res_aa_three)
                save_options()

            elif (aaType == "2" and (codonDetection == "n" or codonDetection == "N" or codonDetection == "no" or codonDetection == "No" or codonDetection == "NO")):
                res_aa_zero = fastaWizard.template_to_aminoacid[0]
                print("\nThe three-letter amino acid sequence of your given sequence is: ")
                print(res_aa_zero)
                save_options()

            elif (aaType == "2" and (codonDetection == "y" or codonDetection == "Y" or codonDetection == "yes" or codonDetection == "Yes" or codonDetection == "YES")):
                res_aa_two = fastaWizard.template_to_aminoacid[2]
                print("\nThe three-letter amino acid sequence of your given sequence is: ")
                print(res_aa_two)
                save_options()

            elif (aaType == "3" and (codonDetection == "n" or codonDetection == "N" or codonDetection == "no" or codonDetection == "No" or codonDetection == "NO")):
                res_aa_both_one = fastaWizard.template_to_aminoacid[1]
                fastaWizard = FastaWizard()
                res_aa_both_zero = fastaWizard.template_to_aminoacid[0]
                print("\nThe amino acid sequences of your given sequence are: ")
                print(res_aa_both_one)
                print(res_aa_both_zero)
                save_options()

            elif (aaType == "3" and (codonDetection == "y" or codonDetection == "Y" or codonDetection == "yes" or codonDetection == "Yes" or codonDetection == "YES")):
                res_aa_both_three = fastaWizard.template_to_aminoacid[3]
                fastaWizard = FastaWizard()
                res_aa_both_two = fastaWizard.template_to_aminoacid[2]
                print("\nThe amino acid sequences of your given sequence are: ")
                print(res_aa_both_three)
                print(res_aa_both_two)
                save_options()

        elif (subprocessOne == "5"):
            res_all_complementary = fastaWizard.template_to_complementary()
            res_all_codon = fastaWizard.template_to_codon()
            res_all_anticodon = fastaWizard.template_to_anticodon()
            fastaWizard = FastaWizard()
            res_all_aa_single = fastaWizard.template_to_aminoacid[1]
            fastaWizard = FastaWizard()
            res_all_aa_triple = fastaWizard.template_to_aminoacid[0]
            print("\nThe complementary of your given sequence is: ")
            print(res_all_complementary)
            print("\nThe codon (mRNA sequence) of your given sequence is: ")
            print(res_all_codon)
            print("\nThe anticodon (tRNA sequence) of your given sequence is: ")
            print(res_all_anticodon)
            print("\nThe single-letter amino acid sequence of your given sequence is: ")
            print(res_all_aa_single)
            print("\nThe three-letter amino acid sequence of your given sequence is: ")
            print(res_all_aa_triple)
            time.sleep(5)
            save_options()

        elif (subprocessOne == "q"):
            continue

        else:
            print("\033[31mERROR: Please choose one option or return to main menu... Do not enter unordinary strings or integers to input box...\033[0m")
            time.sleep(4)
            pass


    elif (process == "2"):
        sequence_options()
        print("""
\033[32m**********************************
1) G-C content of sequence
        
2) A-T content of sequence
        
3) A, T, G, C ratios of DNA sequence

4) A, U, G, C ratios of RNA sequence

5) Full content analysis
        
Please press 'q' to return main menu
**********************************\033[0m\n\n
        """)
        time.sleep(1)
        subprocessTwo = input("Please select to subprocess you would like to do: ")

        if (subprocessTwo == "1"):
            res_gc_content = fastaWizard.gc_content
            print("G-C concentration is " + str(res_gc_content) + "% of your given sequence")
            save_options()

        elif (subprocessTwo == "2"):
            res_at_content = fastaWizard.at_content
            print("A-T concentration is " + str(res_at_content) + "% of your given sequence")
            save_options()

        elif (subprocessTwo == "3"):
            res_dna_aratio = fastaWizard.adenine_ratio
            res_dna_tratio = fastaWizard.thymine_ratio
            res_dna_gratio = fastaWizard.guanine_ratio
            res_dna_cratio = fastaWizard.cytosine_ratio
            print("Adenine ratio = " + "%" + str(res_dna_aratio))
            print("Thymine ratio = " + "%" + str(res_dna_tratio))
            print("Guanine ratio = " + "%" + str(res_dna_gratio))
            print("Cytosine ratio = " + "%" + str(res_dna_cratio))
            save_options()

        elif (subprocessTwo == "4"):
            res_rna_aratio = fastaWizard.adenine_ratio
            res_rna_uratio = fastaWizard.uracile_ratio
            res_rna_gratio = fastaWizard.guanine_ratio
            res_rna_cratio = fastaWizard.cytosine_ratio
            print("Adenine ratio = " + "%" + str(res_rna_aratio))
            print("Uracile ratio = " + "%" + str(res_rna_uratio))
            print("Guanine ratio = " + "%" + str(res_rna_gratio))
            print("Cytosine ratio = " + "%" + str(res_rna_cratio))
            save_options()

        elif (subprocessTwo == "5"):
            res_all_aratio = fastaWizard.adenine_ratio
            res_all_tratio = fastaWizard.thymine_ratio
            res_all_gratio = fastaWizard.guanine_ratio
            res_all_cratio = fastaWizard.cytosine_ratio
            res_all_gc = fastaWizard.gc_content
            res_all_at = fastaWizard.at_content
            print("Adenine ratio = " + "%" + str(res_all_aratio))
            print("Thymine ratio = " + "%" + str(res_all_tratio))
            print("Guanine ratio = " + "%" + str(res_all_gratio))
            print("Cytosine ratio = " + "%" + str(res_all_cratio))
            print("\nG-C concentration is " + str(res_all_gc) + "% of your given sequence")
            print("\nA-T concentration is " + str(res_all_at) + "% of your given sequence")
            save_options()

        elif (subprocessTwo == "q"):
            continue

        else:
            print("\033[31mERROR: Please choose one option or return to main menu... Do not enter unordinary strings or integers to input box...\033[0m")
            time.sleep(4)
            pass


    elif(process == "3"):
        subprocessOne = ""
        subprocessTwo = ""
        sequence_options()

        for rotate in range(101):
            sys.stdout.write("\r Similarity calculating ... %{}".format(rotate))
            sys.stdout.flush()
            time.sleep(0.1)

        similarity = fastaWizard.sequence_similarity()

        if (0 <= similarity <= 20):
            print("""
\033[32moooo\033[0m\033[31moooooooooooooooo\033[0m
\033[32m||||\033[0m\033[31m||||||||||||||||\033[0m
\033[32moooo\033[0m\033[31moooooooooooooooo\033[0m
            """)

        elif (20 < similarity <= 40):
            print("""
\033[32moooooooo\033[0m\033[31moooooooooooo\033[0m
\033[32m||||||||\033[0m\033[31m||||||||||||\033[0m
\033[32moooooooo\033[0m\033[31moooooooooooo\033[0m
            """)

        elif (40 < similarity <= 60):
            print("""
\033[32moooooooooooo\033[0m\033[31moooooooo\033[0m
\033[32m||||||||||||\033[0m\033[31m||||||||\033[0m
\033[32moooooooooooo\033[0m\033[31moooooooo\033[0m
            """)

        elif (60 < similarity <= 80):
            print("""
\033[32moooooooooooooooo\033[0m\033[31moooo\033[0m
\033[32m||||||||||||||||\033[0m\033[31m||||\033[0m
\033[32moooooooooooooooo\033[0m\033[31moooo\033[0m
            """)

        elif (80 < similarity <= 99):
            print("""
\033[32mooooooooooooooooooo\033[0m\033[31mo\033[0m
\033[32m|||||||||||||||||||\033[0m\033[31m|\033[0m
\033[32mooooooooooooooooooo\033[0m\033[31mo\033[0m
            """)

        elif (99 < similarity):
            print("""
\033[32moooooooooooooooooooo\033[0m
\033[32m||||||||||||||||||||\033[0m
\033[32moooooooooooooooooooo\033[0m
             """)

        print("\nThe nucleotide similarity between your sequences are " + str(similarity) + " %")
        save_options()


    elif (process == "q"):
        print("Exiting the program...")
        time.sleep(4)
        break


    else:
        print("\033[31mERROR: Please choose one option or quit program... Do not enter unordinary strings or integers to input box...\033[0m")
        time.sleep(4)
        continue