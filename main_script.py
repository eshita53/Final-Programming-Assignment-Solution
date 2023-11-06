'''
This is the main module where all the class  and analysis calling will be done
'''
import argparse
from matplotlib import pyplot as plt
import yaml
import expression_data_class as EDC
import statistica_analysis as SA
import custom_exceptions as E


def get_config():
    '''
    This function calls the config.yaml file
    and fetches the configuration files
    '''
    with open('config.yaml', 'r', encoding='UTF-8') as stream:
        config = yaml.safe_load(stream)
    return config


def args_parser():
    '''
    This is the args parser function, described
    which parameters are required to execute
    the script
    '''
    parser = argparse.ArgumentParser(
        description=" This takes the gene names ")
    parser.add_argument('-gn', '--gene_name', type=str,
                        help='Name of a gene')

    return parser


def scater_plot(x, y):
    '''
    This function takes the p-values and idfferential expression for all genes
    list where x cordinates is the differentials_expression_all_genes
    and P_values_all_genes is the y cordinates
    '''
    plt.scatter(x, y)
    plt.xlabel("Differential gene Expressions Normal versus HCC")
    plt.ylabel("-log adjusted P_Value")
    plt.show()


def main():
    '''
    This is the main module 
    '''
    config = get_config()
    file_name = config['file_name']
    cancer_state = config['cancer_state']
    expression_data_class = EDC.ExpressionDataClass(file_name, cancer_state)
    statistical_analysis_class = SA.StatisticalAnalysis(expression_data_class)

    try:
        args = args_parser().parse_args()

        # taking the gene name from args parser
        gene_name = args.gene_name

        # calling the get_gene_names function from the expression data class object
        gene_names_list = expression_data_class.get_gene_names()

        # calling the get_expression function from the expression data class object using
        # the gene name from args parser and
        gene_expression_values = expression_data_class.get_expression(
            gene_name)

        # calling the calculate_mean_expression function from statistical analysis object
        # to get a mean expression value across all samples

        mean_expression_for_given_gene = statistical_analysis_class.calculate_mean_expression(
            gene_name)

        # calling the calculate_p_value function from statistical analysis object
        # to get a p value of  all genes as a dictionary where gene_name is the key
        # p_value is the value.
        gene_name_p_values_dict = statistical_analysis_class.calculating_the_p_value()

        # calling differential_expression_for_all_genes functions to get
        # the differentials expresions for all genes
        differentials_expression_all_genes_dictionary = statistical_analysis_class.\
            differential_expression_for_all_genes()

        # from statistical_analysis class object, calling the most expressed genes
        # fuctions to get the expressed gene names and taking the top 10
        # by using list comprehensions and generator expressions
        most_expressed_genes = statistical_analysis_class.\
            most_expressed_genes_in_normal_vs_cancerous_state()
        top_10_most_expressed_genes_normal_vs_cancerous = [
            next(most_expressed_genes) for i in range(10)]

        # from statistical_analysis class object, calling the less expressed genes
        # fuctions to get the expressed gene names and taking the top 10
        # by using list comprehensions and generator expressions

        less_expressed_genes = statistical_analysis_class.\
            less_expressed_genes_in_normal_vs_cancerous_state()
        top_10_less_expressed_genes_normal_vs_cancerous = [
            next(less_expressed_genes) for i in range(10)]

        # calling the scatter_lpot function to plot
        # the graph where p values for all genes will be in the y axis
        # and differentials expressions for all genes will be in the x axis
        y = gene_name_p_values_dict.values()
        x = differentials_expression_all_genes_dictionary.values()
        scater_plot(x, y)

        # writing the output of gene_names_list and  gene_expression_values in the output.txt file
        with open('output.txt', 'w', encoding='UTF-8') as file:
            file.writelines(f'Gene Name list:\n{gene_names_list}\n')
            file.writelines(
                f'Gene expression Values across samples:\n{gene_expression_values}\n')
            file.writelines(
                f'Mean expression for {gene_name} across sample:\
                {mean_expression_for_given_gene}\n')
            file.writelines(
                f'Top 10 most expressed genes in normal versus cance state \
                \n{top_10_most_expressed_genes_normal_vs_cancerous}\n')
            file.writelines(
                f'Top 10 less expressed genes in normal versus cance state \
               \n {top_10_less_expressed_genes_normal_vs_cancerous}\n')
    except E.InvalidInput:
        print('Please insert a valid gene names')
    except Exception as error:
        raise error

    return 0


if __name__ == '__main__':
    main()
