'''
This script does the statistical analysis 
of the expression data class.
'''

import math
from statistics import mean
from scipy import stats


class StatisticalAnalysis:
    '''
    This class do all the statistical analysis of ExpressionDataClass.

    Attribute:
    gene_expression_data: It is the only required attribute of this class. It is the
                          object of Expression Data class
    ################################################################################
                          
    Methods:

    calculate_mean_expression:
    --------------------------
    Parameter:
    gene_name (string): gene name is the only parameters
    Returns (Float):
    Calculates the means expressions values of that particular
    genes across all samples.

    calculate_differential_expression:
    ----------------------------------
    Parameter:
    hcc_state_genes_values(list): List of gene values for hcc state for a single gene
    normal_state_genes_values(list): List of gene values for normal state for a single gene
    Returns(float):
    Returns the logarithem2 based of dividations of mean values of hcc and normal state. It's a float number

    differential_expression_for_all_genes:
    -------------------------------------
    Parameter:
    It does not take any parameter. Use the class attribute
    Returns(Dictionary):
    It returns a dictionary where key is gene name and value is the
    differentials gene expression

    most_expressed_genes_in_normal_vs_cancerous_state:
    --------------------------------------------------
    Parameter:
    It does not take any parameter. Use the class attribute
    Yields(string):
    Each time it yileds a gene name

    less_expressed_genes_in_normal_vs_cancerous_state:
    --------------------------------------------------
    Parameter:
    It does not take any parameter. Use the class attribute
    Yields(string):
    Each time it yileds a gene name

    calculating_the_p_value:
    -----------------------
    Parameter:
    It does not take any parameter. Use the class attribute
    Returns(dictionary):
    It returns a dictionary where key is gene name and value is the
    p_value

    '''

    def __init__(self, gene_expression_data):
        '''
        This methods is called when  StatisticalAnalysis object is created

        attribute: 
        gene_expression_data: It is the instance of the ExpressionDataClass and is the
                              required attribute
        '''

        self.gene_expression_data = gene_expression_data

    def calculate_mean_expression(self, gene_name):
        '''
        This method takes the gene name. Get the gene expression data 
        across all samples from gene expression object instances. and returns the 
        mean expression values.
        '''
        # get the expression values from the gene expression objects
        expression_values = self.gene_expression_data.get_expression(gene_name)
        return mean(expression_values)

    def calculate_differential_expression(self, hcc_state_expression_values,
                                          normal_state_expression_values):
        '''
        This method takes the list of HCC state and normal state gene
        expression values. Calculate the log2 based differentials expressions between them
        and yields the result

        '''
        hcc_state_mean_expression_values = mean(hcc_state_expression_values)
        normal_state_mean_expression_values = mean(
            normal_state_expression_values)

       # used log2-based formula to get the differatials gene expressions

        return math.log2(hcc_state_mean_expression_values / normal_state_mean_expression_values)

    def differential_expression_for_all_genes(self):
        '''
        This methods calculate the differential expression for all the genes
        and store it into a dictionary where gene name is the key
        and differentials expression is the value. And retruns that dictionary
        '''
        gene_name_with_differential_expression_dictionary = {}
        for gene_name, liver_class_objects in \
                self.gene_expression_data.prepared_data_dictionary.items():

            # get the hcc state gene values list from the liver data objects object
            # using list comprehensions where
            # liver data object hcc gene values is not zero.

            hcc_state_gene_values_list = [
                liver_data_object.hcc_state_gene_value for liver_data_object
                in liver_class_objects if liver_data_object.hcc_state_gene_value != 0]

            # get the normal state gene values list from the liver data objects object
            # using list comprehensions where
            # liver data object normal gene values is not zero.

            normal_state_gene_values_list = [
                liver_data_object.normal_state_gene_value for liver_data_object
                in liver_class_objects if liver_data_object.normal_state_gene_value != 0]

            differential_expression = self.calculate_differential_expression(
                hcc_state_gene_values_list, normal_state_gene_values_list)
            # p_value = self.calculating_the_p_value(hcc_state_gene_values_list, normal_state_gene_values_list)
            gene_name_with_differential_expression_dictionary[gene_name] = differential_expression

        return gene_name_with_differential_expression_dictionary

    def most_expressed_genes_in_normal_vs_cancerous_state(self):
        '''
        This methods calles the differential_expression_for_all_genes
        functions to get the dictionary where gene_names is the key
        and differentials equations is the values. The dictionary is sorted according to the 
        reverse orderr to get the most expressed genes. Used yield to send the keys it requires
        to show the most expressed gene names.

        '''

        genes_with_differential_expressions = self.differential_expression_for_all_genes()

        # Sorting the dictionary in decending order and keeping it into a dictionary

        most_expressed_genes = dict(
            sorted(genes_with_differential_expressions.items(), key=lambda e: e[1], reverse=True))

        # Yielding the key when the methods is called
        for keys in most_expressed_genes.keys():
            yield keys

    def less_expressed_genes_in_normal_vs_cancerous_state(self):
        '''
        This methods calles the differential_expression_for_all_genes
        functions to get the dictionary where gene_names is the key
        and differentials equations is the values. The dictionary is sorted according to the 
        ascending order to get the less expressed genes. Used yield to send the keys it requires
        to show the less expressed gene names.

        '''

        genes_with_differential_expressions = self.differential_expression_for_all_genes()

        # Sorting the dictionary and keeping it into dictionary
        less_expressed_genes = dict(
            sorted(genes_with_differential_expressions.items(), key=lambda e: e[1])) 
        
        # Yielding the key when the methods is called
        for keys in less_expressed_genes.keys():
            yield keys

    def calculating_the_p_value(self):
        '''This methods calculate the p_value assuming the
        data has been normally distributed, and returns dictionary with
        gene names as a key and p_value as a value.
        '''

        gene_name_t_test_p_value_dict = {}
        for gene_name, liver_class_objects in \
                self.gene_expression_data.prepared_data_dictionary.items():

            # get the hcc state gene values list from the liver data objects object
            # using list comprehensions where
            # liver data object hcc gene values is not zero.

            hcc_state_gene_values_list = [
                liver_data_object.hcc_state_gene_value for liver_data_object
                in liver_class_objects if liver_data_object.hcc_state_gene_value != 0]

            # get the normal state gene values list from the liver data objects object
            # using list comprehensions where
            # liver data object normal gene values is not zero.

            normal_state_gene_values_list = [
                liver_data_object.normal_state_gene_value for liver_data_object
                in liver_class_objects if liver_data_object.normal_state_gene_value != 0]
            _, p_value = stats.ttest_ind(
                hcc_state_gene_values_list, normal_state_gene_values_list, equal_var=True)

          # Here did the minus natural logarithm to get a nice graph as used log2 based
            # in calculating differentials equations

            gene_name_t_test_p_value_dict[gene_name] = -math.log(p_value)
        return gene_name_t_test_p_value_dict
