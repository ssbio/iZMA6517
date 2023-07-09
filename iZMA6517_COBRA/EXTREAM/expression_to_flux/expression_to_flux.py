import pandas
import xlrd
import re

#import data, update the input file name that contains gene expression data
input_workbook1 = xlrd.open_workbook('data.xlsx')

#loading data into pandas dataframe for fetching corresponding gene expression value
df = pandas.read_excel('data.xlsx', sheet_name=0)
df = df.set_index('Genes')

#output1=open('rxn_expression_.txt','w') # change the output file name for each column in the input excel file


# first sheet contains the expression data, second sheet contains the GPR
input_worksheet1 = input_workbook1.sheet_by_index(0) # 0 is the index for the first sheet
number_of_rows1 = input_worksheet1.nrows
number_of_columns1 = input_worksheet1.ncols

# Sheet containing GPR info
input_worksheet2 = input_workbook1.sheet_by_index(1) # 1 is the index for the second sheet
number_of_rows2 = input_worksheet2.nrows
number_of_columns2 = input_worksheet2.ncols


exp_value_column = 'Normalized_expression'
output_file_name = 'rxn_expression_values.txt'
output = open(output_file_name,'w')
exception_list = open("rxn_expression_exception.txt", 'w')


def calculate_value(item):
    print(item)
    items = item.group(0)
    print(f"Parentheses Items: {items}")
    items = items.strip()
    items = items.replace('(','').replace(')','')
    item_list = []
    value_list = []
    if 'and' in items:
        items = items.split('and')
        # print(items)
        for item in items:
            item = item.strip()
            item_list.append(item)
            try:
                value = df.at[item, exp_value_column]
                print(f"Value Found. {item} -----> {value}")
                value_list.append(value)
            except KeyError:
                value_list.append(0.00)

        print(f"Final Item List: {item_list}")        
        print(f"Associated Value List: {value_list}")
        print(f"Returned Value: {min(value_list)}")
        return str(min(value_list))

    elif 'or' in items:
        items = items.split('or')
        print(items)
        for item in items:
            item = item.strip()
            item_list.append(item)
            try:
                value = df.at[item, exp_value_column]
                print(f"Value Found. {item} -----> {value}")
                value_list.append(value)
            except KeyError:
                value_list.append(0.00)

        print(f"Final Item List: {item_list}")        
        print(f"Associated Value List: {value_list}")
        print(f"Returned Value: {sum(value_list)}")
        return str(sum(value_list))
    else:
        return 'YOU DECIDE!!!!!'


for row in range(1, number_of_rows2):

    #reaction 
    rxn = str((input_worksheet2.cell(row, 0).value))

    #gpr
    gpr = str((input_worksheet2.cell(row, 1).value))

    #reaction with NO condition
    if 'and' not in gpr and 'or' not in gpr:
        try:
            value = df.at[gpr, exp_value_column]
        except KeyError:
            value = 1000
        solo_gpr = "'" + rxn + "'" + '\t' + str(value)
        output.write(solo_gpr)
        output.write('\n')

    #reaaction with ONLY AND condition
    elif 'and' in gpr and 'or' not in gpr:
        gpr = gpr.replace('(','').replace(')','')
        and_split_items = gpr.split('and')
        and_split_values = []
        for item in and_split_items:
            item = item.strip()
            try:
                value = df.at[item, exp_value_column]
            except KeyError:
                value = 0
            and_split_values.append(value)

        and_gpr ="'" + rxn + "'" + '\t' + str(min(and_split_values))
        output.write(and_gpr)
        output.write('\n')

    #reaaction with ONLY OR condition
    elif 'or' in gpr and 'and' not in gpr:
        gpr = gpr.replace('(','').replace(')','')
        or_split_items = gpr.split('or')
        or_split_values = []
        for item in or_split_items:
            item = item.strip()
            try:
                value = df.at[item, exp_value_column]
            except KeyError:
                value = 0
            or_split_values.append(value)

        or_gpr ="'" + rxn + "'" + '\t' + str(sum(or_split_values))
        output.write(or_gpr)
        output.write('\n')


    elif 'and' in gpr and 'or' in gpr:
        print("*****************************************************************************")
        print(f"Initial GPR: {gpr}")
        
        gpr_items = re.sub("\(.*?\)", calculate_value, gpr)
        gpr_items = gpr_items.strip()
        print(f"Replaced Equation: {gpr_items}")
        

        
        def calculate_final_value(gpr_items):
            cleaned_gpr_items = []
            gpr_values = []
            for item in gpr_items:
                item = item.strip()
                cleaned_gpr_items.append(item)
            for item in cleaned_gpr_items:
                if item[0].isdigit():
                    print(f"{item} -----> Number")
                    try:
                        gpr_values.append(float(item))
                    except ValueError:
                        gpr_values.append(0.00)
                        exception_value = "'" + rxn + "'" + '\t' + gpr
                        exception_list.write(exception_value)
                        exception_list.write('\n')
                        
                else:
                    print(f"{item} -----> GPR")
                    try:
                        value = df.at[item, exp_value_column]
                        print(f"Value Found. {item} -----> {value}")
                        gpr_values.append(float(value))
                    except KeyError:
                        gpr_values.append(float(0.00))
                        
            print(cleaned_gpr_items)
            print(gpr_values)
            return gpr_values

        final_value = float()

        if 'and' in gpr_items:
            gpr_items = gpr_items.split('and')
            print(gpr_items)
            final_value = min(calculate_final_value(gpr_items))
            print(f"Final Value: {final_value}")
        elif 'or' in gpr_items:
            gpr_items = gpr_items.split('or')
            print(gpr_items)
            final_value = sum(calculate_final_value(gpr_items))
            print(f"Final Value: {final_value}")


        complex_gpr = "'" + rxn + "'" + '\t' + str(final_value)
        print(complex_gpr)
        output.write(complex_gpr)
        output.write('\n')
