import numpy as np
import random
import tkinter as tk

class interaction_approximator():
    """

    This class is designed to hold all functions designed to develop predictions of ssDNA aptamer sequences.

    """

    def __init__(self):

        """
        One dictionary hold data pertaining to bond probability and hydrogen bonding.
        Each amino acid is referenced by one-letter code, existing as keys.
        As values, data pertaining to each nucleotide in relation to the given amino.

        This dictionary holds the probability data in respect to bonds which may occur:
        the first value is the one letter nucleotide code used to develop a predicted string,
        the second element is the pseudopairing probability and quantity,
        the third element is the hydrogen bond probability and quantity.

        """

        self.peptide_collection_text_box = ''

        # the value is in units of kcal/mol
        self.hydrogen_bond_strength = 0.0000000001

        # Each dictionary holds the probability data in respect to hydrogen bonds which may occur
        # the first value is the one letter nucleotide code used to develop a predicted string
        # the second element is the pseudopairing probability and quantity
        # the third element is the hydrogen bond probability and quantity

        # this dictionary is a dictionary of dictionaries -- the first key is a one letter amino
        # acid code -- the second key is a nucleotide which allows access to probability associated values
        self.association_probability_amino_to_nucleotide_dict_dict = {
            'S': {'A': ((1, 0.4), (6, 2.3)), 'G': ((0, 0), (1, 0.4)), 'C': ((0, 0), (3, 12.5))},
            'N': {'A': ((10, 3.8), (4, 1.5)), 'G': ((0, 0), (1, 1.2)), 'C': ((1, 4.2), (0, 0)), 'T': ((2, 8.3), (2, 8.3))},
            'Q': {'A': ((2, 0.8), (4, 1.5)), 'G': ((0, 0), (2, 2.4)), 'T': ((4, 16.7), (1, 4.2))},
            'D': {'A': ((7, 2.6), (10, 3.8)), 'G': ((28, 34.1), (2, 2.4)), 'T': ((0, 0), (1, 4.2))},
            'E': {'A': ((1, 0.4), (3, 1.1)), 'G': ((2, 2.4), (7, 8.5)), 'T': ((0, 0), (1, 4.2))},
            'K': {'T': ((0, 0), (1, 4.2))},
            'T': {'C': ((1, 4.2), (0, 0))},
            'R': {'C': ((3, 12.5), (0, 0))}
        }

        # this dictionary holds data related to each nucleotide
        self.association_probability_nucleotide_to_amino_dict_dict = {
            'A': {'S': ((1, 0.4), (6, 2.3)), 'N': ((10, 3.8), (4, 1.5)), 'Q': ((2, 0.8), (4, 1.5)), 'D': ((7, 2.6), (10, 3.8)),
                  'E': ((1, 0.4), (3, 1.1))},
            'G': {'S': ((0, 0), (1, 0.4)), 'N': ((0, 0), (1, 1.2)), 'Q': ((0, 0), (2, 2.4)), 'D': ((28, 34.1), (2, 2.4)), 'E': ((2, 2.4), (7, 8.5))},
            'C': {'S': ((0, 0), (3, 12.5)), 'N': ((1, 4.2), (0, 0)), 'T': ((1, 4.2), (0, 0)), 'R': ((3, 12.5), (0, 0))},
            'T': {'N': ((2, 8.3), (2, 8.3)), 'Q': ((4, 16.7), (1, 4.2)), 'D': ((0, 0), (1, 4.2)), 'E': ((0, 0), (1, 4.2)), 'K': ((0, 0), (1, 4.2))}
        }

        # the primer sequences must be scored to aid in reducing false positives and obfuscating results.
        # will be used to develop a minimum score to test sequences against for study validity
        self.null_pseudopair_strength_and_peptide_fragment = (0, '')
        self.p5 = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
        self.p7 = "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"
        self.null_model_aptamer_sequence = self.p5 + self.p7
        self.null_pseudopair_bond_strength = 0
        self.null_pseudopair_bond_probability = 0

        self.target_peptide = "MSKGEELFTG VVPILVELDG DVNGHKFSVS GEGEGDATYG KLTLKFICTT GKLPVPWPTL VTTFSYGVQC FSRYPDHMKQ HDFFKSAMPE GYVQERTIFF KDDGNYKTRA EVKFEGDTLV NRIELKGIDF KEDGNILGHK LEYNYNSHNV YIMADKQKNG IKVNFKIRHN IEDGSVQLAD HYQQNTPIGD GPVLLPDNHY LSTQSALSKD PNEKRDHMVL LEFVTAAGIT HGMDELYK"
        self.peptide_fragment_list = []

        self.peptide_fragment_and_corresponding_aptamer_tup_list = []
        self.desired_aptamer_length = 50

        self.threshold_value = 0

        # set to false to alter sort contingent variable
        self.sort_by_bond_strength = True

        self.programatc_iterations = 10

    def peptide_parser(self, desired_length):
        """ This function is designed to parse the given peptide into all possible sequences of the desired length.
        These fragments will be considered as possible locations for aptamer interaction.


        Input: peptide sequence
        Output: list of possible interaction locations
        """

        target_peptide_stripped = self.target_peptide.replace(" ", '')

        peptide_fragment_list = []
        for index in range(0, len(target_peptide_stripped) - desired_length):
            peptide_fragment_list.append(target_peptide_stripped[index:index + desired_length])

        return peptide_fragment_list

    def null_model_aptamer_value_calculator(self):
        """The null model in the context of this algorithm is the result of prediction using the
        PCR primers alone, represented by the P5 and P7 primers. As such the null model is deterministic.

        Values are formed by summing the number of bonds and probabilities associated with both hydrogen bonding
        and pseudopairing. This develops a total bonding score.

        The maximum of the values developed by evaluation of all fragments is the null value. The sequence which
        develops the maximum total bonding score is identified and saved.

        Input: The primer sequence.
        Output: A value representing the score of the primers alone.
        """

        null_model_peptide_fragment_list = self.peptide_parser(len(self.p5) + len(self.p7))
        null_model_aptamer = self.p5 + self.p7

        null_model_sum_list = []
        for peptide_fragment in null_model_peptide_fragment_list:
            null_model_score_list = []
            for i in range(len(peptide_fragment)):
                if peptide_fragment[i] in self.association_probability_nucleotide_to_amino_dict_dict[null_model_aptamer[i]].keys():
                    # condense values
                    null_model_pseudopair_data = \
                    self.association_probability_nucleotide_to_amino_dict_dict[null_model_aptamer[i]][peptide_fragment[i]][0]
                    null_model_hydrogen_bond_data = \
                    self.association_probability_nucleotide_to_amino_dict_dict[null_model_aptamer[i]][peptide_fragment[i]][1]

                    bond_sum_data = (null_model_pseudopair_data[0] + null_model_hydrogen_bond_data[0],
                                     null_model_pseudopair_data[1] + null_model_hydrogen_bond_data[1])

                    null_model_score_list.append(bond_sum_data)

            null_model_sum_tuple = sum(i[0] for i in null_model_score_list), sum(i[1] for i in null_model_score_list)
            null_model_sum_list.append((self.p5 + self.p7, null_model_sum_tuple, peptide_fragment))

        # deterministic null model value used as threhold to determine viability of aptamer candidates by sum bond probability
        self.threshold_value = sum(j[1][1] for j in null_model_sum_list) / len(null_model_sum_list)

    def theoretical_aptamer_predictor(self):
        """
        In developing an aptamer from the peptide sequence, the aptamer sequnce is treated as a hidden path.
        The peptide sequence is treated as the empirical path.

        As the peptide path is all that is known, the peptide to amino dictionary dictionary will be used.

        Nucleotides will be chosen by use of the probability values as weights for the random choice library function.

        As is similar to the means of developing the null model data, a total bonding score will be developed for each possible interaction.

        Input: peptide sequence.
        Output: aptamer sequences paired with total bonding scores.

        """

        peptide_fragment_list = self.peptide_parser(self.desired_aptamer_length)

        for peptide_fragment in peptide_fragment_list:
            theoretical_aptamer = ''
            bond_prob_sum = 0
            total_prob_bonds = 0

            for i in range(0, len(peptide_fragment)):
                nucleotide_associated_probs_sum_list = []
                possible_base_list = []
                sum_tup_list = []

                if peptide_fragment[i] in self.association_probability_amino_to_nucleotide_dict_dict.keys():
                    for possible_nucleotide in self.association_probability_amino_to_nucleotide_dict_dict[peptide_fragment[i]]:
                        # extract sum probabilities into tuple
                        prob_bond_tup_list = self.association_probability_amino_to_nucleotide_dict_dict[peptide_fragment[i]][possible_nucleotide]

                        probability_sum_tup = possible_nucleotide, sum(j[0] for j in prob_bond_tup_list), sum(j[1] for j in prob_bond_tup_list)

                        sum_tup_list.append(probability_sum_tup)

                        possible_base_list.append(probability_sum_tup[0])

                        nucleotide_associated_probs_sum_list.append(probability_sum_tup)

                    base_probability_weight_tup = tuple([i[1] for i in nucleotide_associated_probs_sum_list])

                    # use probability sums as weights to select positionally corresponding bases
                    selected_base = random.choices(possible_base_list, weights=base_probability_weight_tup, k=1)[0]
                    for sum_tup in sum_tup_list:
                        if sum_tup[0] == selected_base:
                            bond_prob_sum += sum_tup[1]
                            total_prob_bonds += sum_tup[2]

                    theoretical_aptamer += selected_base

                else:

                    theoretical_aptamer += '-'

            theoretical_aptamer_tup = (theoretical_aptamer, peptide_fragment, bond_prob_sum, total_prob_bonds)

            self.threshhold_checker(theoretical_aptamer_tup)

    def threshhold_checker(self, theoretical_aptamer_tup):

        if theoretical_aptamer_tup[2] > self.threshold_value:
            self.peptide_fragment_and_corresponding_aptamer_tup_list.append((theoretical_aptamer_tup[0], theoretical_aptamer_tup[1],
                                                                             theoretical_aptamer_tup[2], round(theoretical_aptamer_tup[3], 3)))

    def print_out_results(self):
        self.peptide_fragment_and_corresponding_aptamer_tup_list.sort(key=lambda x: x[2], reverse=True)
        for tup in self.peptide_fragment_and_corresponding_aptamer_tup_list:
            print(tup)

    def driver(self):
        """The purpose of this function is to drive iterations of the program.
        """
        print("programatic iterations: ", self.programatc_iterations)
        print(" ")

        # launch input collection/results GUI

        # find bonding strength and probability scores for primer sequences p5 and p7
        self.null_model_aptamer_value_calculator()

        for i in range(0, self.programatc_iterations):
            self.theoretical_aptamer_predictor()

        self.print_out_results()

    def peptide_text_box_collect(self):

        # Insert
        self.target_peptide = self.peptide_collection_text_box.get()
        print(self.target_peptide)

    def GUI(self):

        root_window = tk.Tk()
        root_window.geometry("1000x750")

        root_window.wm_title("Aptamer Predictor by Quin Lamothe for Bernick Lab UCSC, 2022")

        iterant_entry = tk.Entry()


        peptide_entry_label = tk.Label(text="Target Peptide Entry:")
        self.peptide_collection_text_box = tk.Entry()
        peptide_entry_label.pack()
        self.peptide_collection_text_box.pack()

        peptide_collection_button = tk.Button(root_window, text="Enter", command=self.peptide_text_box_collect)
        peptide_collection_button.pack()

        driver_button = tk.Button(root_window, text="Process", command=self.driver)
        driver_button.pack()
        root_window.mainloop()


def main():

    class_access = interaction_approximator()

    class_access.GUI()

if __name__ == '__main__':
    main()