'''
This script represnts the Liver
Data Class.
'''


class LiverData:
    '''
    This class represent the data for a
    gene.

    Attribute:
    sample_name: Represents the sample name
    hcc_state_gene_value: This represents the cancerous state gene value
    norrmal_state_gene_value: Represents the gene value for the normal state

    '''

    def __init__(self, sample_name, hcc_state_gene_value, normal_state_gene_value):
        '''
        This method will be called when an object is created. 
        And it require three attributes for a successfull
        class object creation
        '''
        self.sample_name = sample_name
        self.hcc_state_gene_value = hcc_state_gene_value
        self.normal_state_gene_value = normal_state_gene_value

    def __str__(self):
        '''
        This methods returns the all class attributes 
        '''
        return f'Sample Name:{self.sample_name} ' \
               f'HCC State Gene Values:{self.hcc_state_gene_value} '\
               f'Normal State Gene Values:{self.normal_state_gene_value}'
