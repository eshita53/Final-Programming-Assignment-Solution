'''
This script represnts the Expression 
Data Class
'''

from collections import defaultdict
import liver_data_class as LV


class ExpressionDataClass:
    '''
    This class does all the data processing
    for the given datasets. 

    Attributes:
    file_name: It takes the file name to create the objects of this class. 
               It is the required attribute of the class.
    cancer_state: It is the name of cancer_state and it's another required attribute of this class

    file_raw_data: This attributes is not required to create the class. But it calls the file_reader 
                   methods with the file name to read the file and returns the raw data
                   which can be accessible across the class through this attribute.

    prepared_data_dictionary: This attributes calls the prepare_data methods which 
                              returns a dictinary where gene name is the key and
                              liver data class objects are the values. Liver data
                              class attributes are sample name, hcc state gene values 
                              and normal state gene values.

    ########################################################################

    Methods:

    file_reader
    -----------
    Parameter:
    File name (string): file_name is the only parameters
    Returns (list(list)):  
    It returns a list of file data

    get_gene_names:
    --------------
    Parameter: 
    It does not require any paprameters. It uses the class attributes
    Returns (list):
    A list of gene names

    prepare_data: 
    ------------
    Parameter:
    It does not require any parameters. It uses the class attributes
    Returns (dictionary):
    It returns a dictionary where gene name is the key and
    liver data class objects are the values. Liver data
    class attributes have sample name, hcc state gene values 
    and normal state gene values.    

    get_expression:
    ---------------
    Parameter (string):
    given_gene_name:
    It takes a gene name as it's parameter.
    Returns (list):
    It returns a list of gene expression values across all samples

    '''

    def __init__(self, file_name, cancer_state):
        '''
        This methods is called when an object is created. And it required file_name and
        cancer_state to create the class instance
        '''
        self.cancer_state = cancer_state
        self.file_raw_data = self.file_reader(file_name)
        self.prepared_data_dictionary = self.prepare_data()

    def file_reader(self, file_name):
        '''
        This methods read the file and return the raw data
        '''

        with open(file_name, 'r', encoding='UTF-8') as file:
            liver_data = file.readlines()

        return liver_data

    def prepare_data(self):
        '''
        This methods prepare the raw data of liver file. It 
        converts the data into a dictionary where gene name is the key
        and created a liver data class as the dictinary values where
        sample name, state type and in that state whats the gene value
        is the attribute of that class

        Arguments: Used Liverdata class object as a dictinary to structure
                  data in a more managable way.

        '''

        liver_data_dictionary = defaultdict(list)

        # Calling the get_gene_names functions to get the
        # the list of gene names which can be used as the key
        # names for the liver_data_dictionary.

        gene_names_list = self.get_gene_names()

        # Taking a empty list to store the tuples
        # of sample name, state type and gene values from the raw data

        data_list = []

        # It is taking raw data from the class attributes.
        # Striping the data to remove any extra lines from the end,
        # Spliting the data by commas, and appending it into
        # the data_list as a tuple of sample, state type and gene values.

        for data in self.file_raw_data[1:]:
            data = data.strip()
            sample, state_type, *gene_values = data.split(',')
            data_list.append((sample, state_type, gene_values))

        # Iterating through the gene names list to get the
        # gene names and it's index to use the gene name as the
        # key of liver data dictionary. Using List comprehensions,
        # from the data list taking index wise gene values, converted 
        # the values into float for further mathmatical operations
        # upon depending on the state type, creating a Liver data
        # class with sample name, hcc_type_gene values,
        # normal_state gene_values. when the state from data_list is
        # HCC the normal state gene values will be zeros and it will be
        # vice versa. List of Liver data class objects will be the
        # values of the dictionary.

        for index, gene_name in enumerate(gene_names_list):

            liver_data_dictionary[gene_name] = [
                LV.LiverData(sample_name, float(gene_values[index]), 0)
                if state_type == self.cancer_state
                else LV.LiverData(sample_name, 0, float(gene_values[index]))
                for sample_name, state_type, gene_values in data_list
            ]

        # clearing the memory as we don't need
        # the tuple list of sample_name, state_type, gene_values
        # anymore. The data has been transfered into the dictionary

        del data_list

        return liver_data_dictionary

    def get_gene_names(self):
        '''This methods returns the list of gene names
        from the data'''
        gene_names_list = [names.strip()
                           for names in self.file_raw_data[0].split(",")[2:]]

        return gene_names_list

    def get_expression(self, given_gene_name):
        '''
        This methods takes the gene name and return a list
        of its expression values across different samples
        '''
        # First get values for the given gene name
        # from the class attribute prepared data dictionary
        # Using list comprehensions from the liver data object
        # we take hcc state gene values and normal state gene value
        # when hcc state value is zero.
        #  At the end, get a list of
        # gene expression values for the given gene across all samples

        gene_expression_across_all = [
            liver_data_object.hcc_state_gene_value
            if liver_data_object.hcc_state_gene_value != 0
            else liver_data_object.normal_state_gene_value
            for liver_data_object in self.prepared_data_dictionary[given_gene_name]]

        return gene_expression_across_all

