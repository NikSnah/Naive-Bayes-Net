from bnetbase import Variable, Factor, BN
import csv
import itertools


def normalize(factor):
    '''
    Normalize the factor such that its values sum to 1.
    Do not modify the input factor.

    :param factor: a Factor object. 
    :return: a new Factor object resulting from normalizing factor.
    '''
    # Create a new factor with the same scope
    new_factor = Factor(factor.name + '_norm', factor.get_scope())

    # Copy the values from the original factor
    new_factor.values = list(factor.values)

    # Sum all the values
    total = sum(new_factor.values)
    #print('1', total)
    if total == 0:
        # Avoid division by zero
        return new_factor
    
    # Normalize the values
    new_factor.values = [v / total for v in new_factor.values]

    # print('2', new_factor.print_table())
    return new_factor


def restrict(factor, variable, value):
    '''
    Restrict a factor by assigning value to variable.
    Do not modify the input factor.

    :param factor: a Factor object.
    :param variable: the variable to restrict.
    :param value: the value to restrict the variable to
    :return: a new Factor object resulting from restricting variable to value.
             This new factor no longer has variable in it.

    '''
    # Create new scope without the variable
    new_scope = factor.get_scope()
    # print('1', new_scope)
    if variable not in new_scope:
        return factor  # Variable not in scope, return original factor
    
    new_scope.remove(variable)
    new_factor = Factor(factor.name + '_restrict_' + variable.name, new_scope)

    # Iterate over all assignments in the new scope
    assignments = list(itertools.product(*[var.domain() for var in new_scope]))
    # print('2', assignments)
    for assignment in assignments:
        # Set variable assignments
        for var, val in zip(new_scope, assignment):
            var.set_assignment(val)
        # Set the restricted variable's assignment
        variable.set_assignment(value)
        # Get the value from the original factor
        val = factor.get_value_at_current_assignments()
        # Add value to the new factor
        new_factor.add_value_at_current_assignment(val)

    # print('3', new_factor.print_table())
    return new_factor



def sum_out(factor, variable):
    '''
    Sum out a variable variable from factor factor.
    Do not modify the input factor.

    :param factor: a Factor object.
    :param variable: the variable to sum out.
    :return: a new Factor object resulting from summing out variable from the factor.
             This new factor no longer has variable in it.
    '''
    # Create new scope without the variable
    new_scope = factor.get_scope()

    if variable not in new_scope:
        return factor  # Variable not in scope, return original factor
    

    new_scope.remove(variable)
    new_factor = Factor(factor.name + '_sumout_' + variable.name, new_scope)
    # Iterate over all assignments in the new scope
    assignments = list(itertools.product(*[var.domain() for var in new_scope]))

    for assignment in assignments:
        # Set variable assignments
        for var, val in zip(new_scope, assignment):
            var.set_assignment(val)
        total = 0
        # Sum over all values of the variable being summed out
        for val in variable.domain():
            variable.set_assignment(val)
            total += factor.get_value_at_current_assignments()
        # Add summed value to the new factor
        new_factor.add_value_at_current_assignment(total)

    return new_factor

def multiply(factor_list):
    '''
    Multiply a list of factors together.
    Do not modify any of the input factors. 

    :param factor_list: a list of Factor objects.
    :return: a new Factor object resulting from multiplying all the factors in factor_list.
    '''
    # Union of all variables in the factors
    all_vars = []
    for factor in factor_list:
        for var in factor.get_scope():
            if var not in all_vars:
                all_vars.append(var)
    
   # print('2', all_vars)

    # Create a new factor with the union scope
    new_factor = Factor('multiply_' + '_'.join([f.name for f in factor_list]), all_vars)


    # Iterate over all assignments in the union scope
    assignments = list(itertools.product(*[var.domain() for var in all_vars]))
    #print('3', assignments)
    for assignment in assignments:
        #print('4', assignment)
        for var, val in zip(all_vars, assignment):
            var.set_assignment(val)
        # Multiply values from each factor
        total = 1
        for factor in factor_list:
            total *= factor.get_value_at_current_assignments()
        #print('5', total)
        #print('6', new_factor.get_scope())
        new_factor.add_value_at_current_assignment(total)

    #print ('8', new_factor.print_table())

    
    return new_factor

def ve(bayes_net, var_query, EvidenceVars):
    '''

    Execute the variable elimination algorithm on the Bayesian network bayes_net
    to compute a distribution over the values of var_query given the 
    evidence provided by EvidenceVars. 

    :param bayes_net: a BN object.
    :param var_query: the query variable. we want to compute a distribution
                     over the values of the query variable.
    :param EvidenceVars: the evidence variables. Each evidence variable has 
                         its evidence set to a value from its domain 
                         using set_evidence.
    :return: a Factor object representing a distribution over the values
             of var_query. that is a list of numbers, one for every value
             in var_query's domain. These numbers sum to 1. The i-th number
             is the probability that var_query is equal to its i-th value given 
             the settings of the evidence variables.

    For example, assume that
        var_query = A with Dom[A] = ['a', 'b', 'c'], 
        EvidenceVars = [B, C], and 
        we have called B.set_evidence(1) and C.set_evidence('c'), 
    then VE would return a list of three numbers, e.g. [0.5, 0.24, 0.26]. 
    These numbers would mean that 
        Pr(A='a'|B=1, C='c') = 0.5, 
        Pr(A='a'|B=1, C='c') = 0.24, and 
        Pr(A='a'|B=1, C='c') = 0.26.

    '''
     # Initialize factors
    factors = bayes_net.factors()
    
    # Apply evidence by restricting factors
    for ev_var in EvidenceVars:
        value = ev_var.get_evidence()
        new_factors = []
        for factor in factors:
            if ev_var in factor.get_scope():
                new_factor = restrict(factor, ev_var, value)
                new_factors.append(new_factor)

                # print('1 new_factor', new_factor.print_table())
            else:
                new_factors.append(factor)
                # print('1 factor', factor.print_table())
        factors = new_factors
    
    # Variables to eliminate
    vars_to_eliminate = set(bayes_net.variables())
    vars_to_eliminate.discard(var_query)
    vars_to_eliminate.difference_update(EvidenceVars)
    vars_to_eliminate = list(vars_to_eliminate)

    # print('2', vars_to_eliminate)
    
    # Eliminate variables
    for var in vars_to_eliminate:
        # Factors involving the variable
        factors_with_var = [f for f in factors if var in f.get_scope()]
        if factors_with_var:
            # Remove these factors from the list
            factors = [f for f in factors if var not in f.get_scope()]
            # Multiply factors
            multiplied = multiply(factors_with_var)
            # Sum out the variable
            summed_out = sum_out(multiplied, var)
            factors.append(summed_out)

    # print('3', factors)
    
    # Multiply remaining factors
    if factors:
        result_factor = multiply(factors)
    else:
        # No factors left
        result_factor = Factor('unit', [var_query])
        for val in var_query.domain():
            var_query.set_assignment(val)
            result_factor.add_value_at_current_assignment(1)

    # print('4', result_factor.print_table())
    
    # Normalize the result
    result_factor = normalize(result_factor)

    # print('5', result_factor.print_table())
    return result_factor


def naive_bayes_model(data_file, variable_domains = {"Work": ['Not Working', 'Government', 'Private', 'Self-emp'], "Education": ['<Gr12', 'HS-Graduate', 'Associate', 'Professional', 'Bachelors', 'Masters', 'Doctorate'], "Occupation": ['Admin', 'Military', 'Manual Labour', 'Office Labour', 'Service', 'Professional'], "MaritalStatus": ['Not-Married', 'Married', 'Separated', 'Widowed'], "Relationship": ['Wife', 'Own-child', 'Husband', 'Not-in-family', 'Other-relative', 'Unmarried'], "Race": ['White', 'Black', 'Asian-Pac-Islander', 'Amer-Indian-Eskimo', 'Other'], "Gender": ['Male', 'Female'], "Country": ['North-America', 'South-America', 'Europe', 'Asia', 'Middle-East', 'Carribean'], "Salary": ['<50K', '>=50K']}, class_var = Variable("Salary", ['<50K', '>=50K'])):
    '''
   NaiveBayesModel returns a BN that is a Naive Bayes model that 
   represents the joint distribution of value assignments to 
   variables in the Adult Dataset from UCI.  Remember a Naive Bayes model
   assumes P(X1, X2,.... XN, Class) can be represented as 
   P(X1|Class)*P(X2|Class)* .... *P(XN|Class)*P(Class).
   When you generated your Bayes bayes_net, assume that the values 
   in the SALARY column of the dataset are the CLASS that we want to predict.
   @return a BN that is a Naive Bayes model and which represents the Adult Dataset. 
    '''

    if variable_domains is None:
        variable_domains = {
            "Work": ['Not Working', 'Government', 'Private', 'Self-emp'],
            "Education": ['<Gr12', 'HS-Graduate', 'Associate', 'Professional', 'Bachelors', 'Masters', 'Doctorate'],
            "Occupation": ['Admin', 'Military', 'Manual Labour', 'Office Labour', 'Service', 'Professional'],
            "MaritalStatus": ['Not-Married', 'Married', 'Separated', 'Widowed'],
            "Relationship": ['Wife', 'Own-child', 'Husband', 'Not-in-family', 'Other-relative', 'Unmarried'],
            "Race": ['White', 'Black', 'Asian-Pac-Islander', 'Amer-Indian-Eskimo', 'Other'],
            "Gender": ['Male', 'Female'],
            "Country": ['North-America', 'South-America', 'Europe', 'Asia', 'Middle-East', 'Carribean'],
            "Salary": ['<50K', '>=50K']
        }
    if class_var is None:
        class_var = Variable("Salary", variable_domains["Salary"])

    # Create variables
    variables = {name: Variable(name, domain) for name, domain in variable_domains.items()}

    # Initialize counts
    counts = {}
    for var in variables.values():
        counts[var.name] = {val: 0 for val in var.domain()}
    counts['Joint'] = {}

    # Read data
    with open(data_file, newline='') as csvfile:
        reader = csv.reader(csvfile)
        headers = next(reader, None)  # Skip header row
        data = [row for row in reader]

    # Total number of instances
    total_instances = len(data)

    # Count frequencies
    for row in data:
        assignment = dict(zip(headers, row))
        salary = assignment['Salary']
        counts['Salary'][salary] += 1
        for var_name in headers[:-1]:  # Exclude 'Salary' column
            val = assignment[var_name]
            counts[var_name][val] += 1
            key = (var_name, val, salary)
            counts['Joint'][key] = counts['Joint'].get(key, 0) + 1

    # Create factors
    factors = []
    # P(Salary)
    salary_var = variables['Salary']
    salary_factor = Factor('P(Salary)', [salary_var])
    # print('1', salary_var.domain())
    for val in salary_var.domain():
        prob = counts['Salary'][val] / total_instances
        salary_var.set_assignment(val)
        salary_factor.add_value_at_current_assignment(prob)
    factors.append(salary_factor)

    # P(X|Salary) for each predictor variable X while excluding the salary variable
    for var_name in headers[:-1]:  
        var = variables[var_name]
        factor = Factor(f'P({var_name}|Salary)', [var, salary_var])
        for val in var.domain():
            for s_val in salary_var.domain():
                key = (var_name, val, s_val)
                joint_count = counts['Joint'].get(key, 0)
                salary_count = counts['Salary'][s_val]
                prob = joint_count / salary_count if salary_count > 0 else 0
                var.set_assignment(val)
                salary_var.set_assignment(s_val)
                factor.add_value_at_current_assignment(prob)
        factors.append(factor)

    bayes_net = BN('Naive Bayes', list(variables.values()), factors)
    return bayes_net

    ### READ IN THE DATA
    #input_data = []
   # with open(data_file, newline='') as csvfile:
    #    reader = csv.reader(csvfile)
     #   headers = next(reader, None) #skip header row
    #    for row in reader:
    #        input_data.append(row)

    ### DOMAIN INFORMATION REFLECTS ORDER OF COLUMNS IN THE DATA SET
    #variable_domains = {
    #"Work": ['Not Working', 'Government', 'Private', 'Self-emp'],
    #"Education": ['<Gr12', 'HS-Graduate', 'Associate', 'Professional', 'Bachelors', 'Masters', 'Doctorate'],
    #"Occupation": ['Admin', 'Military', 'Manual Labour', 'Office Labour', 'Service', 'Professional'],
    #"MaritalStatus": ['Not-Married', 'Married', 'Separated', 'Widowed'],
    #"Relationship": ['Wife', 'Own-child', 'Husband', 'Not-in-family', 'Other-relative', 'Unmarried'],
    #"Race": ['White', 'Black', 'Asian-Pac-Islander', 'Amer-Indian-Eskimo', 'Other'],
    #"Gender": ['Male', 'Female'],
    #"Country": ['North-America', 'South-America', 'Europe', 'Asia', 'Middle-East', 'Carribean'],
    #"Salary": ['<50K', '>=50K']
    #}


def compute_probability(bayes_net, e1):
    # variable to query 
    sal_var = bayes_net.get_variable("Salary")

    # set evidence for e1
    for ev in e1:
        ev[0].set_evidence(ev[1])

    # calculate the probability from the bayes_net and variable elimination
    returned_factor = ve(bayes_net, sal_var, [ev[0] for ev in e1])

    return returned_factor.get_value([">=50K"])
    


def explore(bayes_net, question):
    '''    Input: bayes_net---a BN object (a Bayes bayes_net)
           question---an integer indicating the question in HW4 to be calculated. Options are:
           1. What percentage of the women in the data set end up with a P(S=">=$50K"|E1) that is strictly greater than P(S=">=$50K"|E2)?
           2. What percentage of the men in the data set end up with a P(S=">=$50K"|E1) that is strictly greater than P(S=">=$50K"|E2)?
           3. What percentage of the women in the data set with P(S=">=$50K"|E1) > 0.5 actually have a salary over $50K?
           4. What percentage of the men in the data set with P(S=">=$50K"|E1) > 0.5 actually have a salary over $50K?
           5. What percentage of the women in the data set are assigned a P(Salary=">=$50K"|E1) > 0.5, overall?
           6. What percentage of the men in the data set are assigned a P(Salary=">=$50K"|E1) > 0.5, overall?
           @return a percentage (between 0 and 100)
    ''' 

    def set_data(data, type):
        if type == 'e1':
            return [(bayes_net.get_variable("Work"), data[headers.index("Work")]), 
                    (bayes_net.get_variable("Education"), data[headers.index("Education")]), 
                    (bayes_net.get_variable("Occupation"), data[headers.index("Occupation")]), 
                    (bayes_net.get_variable("Relationship"), data[headers.index("Relationship")])]
        else: 
            return  [(bayes_net.get_variable("Work"), data[headers.index("Work")]), 
                    (bayes_net.get_variable("Education"), data[headers.index("Education")]), 
                    (bayes_net.get_variable("Occupation"), data[headers.index("Occupation")]), 
                    (bayes_net.get_variable("Relationship"), data[headers.index("Relationship")]),
                    (bayes_net.get_variable("Gender"), data[headers.index("Gender")])]
                

# Read in the data
    input_data = []
    with open('data/adult-test.csv', newline='') as csvfile:
        reader = csv.reader(csvfile)
        headers = next(reader, None) 
        for row in reader:
            input_data.append(row)

    # calc question answers 
    if question in [1, 2]:
        e1_c = 0
        e2_c = 0
         
        # go through the data and calculate the probability
        for data in input_data:
            e1 = set_data(data, 'e1')
            e2 = set_data(data, 'e2')
            
            p1 = compute_probability(bayes_net, e1)
            p2 = compute_probability(bayes_net, e2)

            #print(p1, p2)
            if question == 1 and data[headers.index("Gender")] == "Female":      
                if p1 > p2:
                    e1_c += 1
                else:
                    e2_c += 1
                    
            elif question == 2 and data[headers.index("Gender")] == "Male":
                if p1 > p2:
                    e1_c += 1
                else:
                    e2_c += 1
            
        #print(p1, p2)
        if question == 1: 
            e2_c = e1_c - e2_c
            return ((e2_c / e1_c) ) * 100
        else:
            return (e1_c / e2_c) * 100

    
    elif question in [3, 4, 5, 6]:
        c = 0
        predict_c = 0

        for data in input_data:
            e1 = set_data(data, 'e1')
            p = compute_probability(bayes_net, e1)
            #print(p)
            if question == 3 and data[headers.index("Gender")] == "Female":                      
                if p > 0.5:
                    predict_c += 1
                    if data[headers.index("Salary")] == ">=50K":
                        c += 1
                #print(c, predict_c)

            elif question == 4 and data[headers.index("Gender")] == "Male":
                if p > 0.5:
                    predict_c += 1
                    if data[headers.index("Salary")] == ">=50K":
                        c += 1
                 #print(c, predict_c)

            elif question == 5 and data[headers.index("Gender")] == "Female":
                if p > 0.5:
                    predict_c += 1
                c += 1
                
            elif question == 6 and data[headers.index("Gender")] == "Male":
                if p > 0.5:
                    predict_c += 1
                c += 1

        #print(c, predict_c)
        if question in [3, 4]:
            return (c / predict_c) * 100
        else:
            return (predict_c / c) * 100

    else: 
        raise ValueError("Invalid question")
    

    
    

if __name__ == '__main__':
    nb = naive_bayes_model('data/adult-train.csv')
    for i in range(1,7):
        print("explore(nb,{}) = {}".format(i, explore(nb, i)))
