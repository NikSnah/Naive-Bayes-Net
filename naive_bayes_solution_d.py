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
    scope = factor.get_scope()
    # Note: normalizing should only occur when the factor has 1 variable. 

    new_name = factor.name + "_normalized"
    new_factor = Factor(new_name, scope)

    domain = list(map(Variable.domain, scope))
    assignments = get_assignments(domain, [])

    total_sum = 0
    factor_values = []
    for assignment in assignments:
        assign_val = factor.get_value(assignment)
        factor_values.append(assign_val)
        total_sum += assign_val


    new_values = []
    for val in factor_values:
        new_values.append(val/total_sum)

    for i in range(len(assignments)):
        assignments[i].append(new_values[i])

    new_factor.add_values(assignments)
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
    scope = factor.get_scope()
    index = scope.index(variable)


    # new scope without variable.
    new_scope = []
    for var in scope:
        if var != variable:
            new_scope.append(var)

    new_name = factor.name + "_restr_" + variable.name + "=" + value

    new_factor = Factor(new_name, new_scope)

    curr_domain = list(map(Variable.domain,scope))
    curr_domain[index] = [value]

    assignments = get_assignments(curr_domain, [], 0)
    for assignment in assignments:
        assignment.append(factor.get_value(assignment))

    for assignment in assignments:
        assignment.pop(index)

    new_factor.add_values(assignments)

    return new_factor
       
        
def get_assignments(domain, curr_assignment=[], index=0):
    assignments = []
    if index == len(domain):
        return [curr_assignment]
    for value in domain[index]:
        new_assignment = list(curr_assignment)
        new_assignment.append(value)
        assignments.extend(get_assignments(domain, new_assignment, index+1))


    return assignments



    

def sum_out(factor, variable):
    '''
    Sum out a variable variable from factor factor.
    Do not modify the input factor.

    :param factor: a Factor object.
    :param variable: the variable to sum out.
    :return: a new Factor object resulting from summing out variable from the factor.
             This new factor no longer has variable in it.
    '''
    scope = factor.get_scope()
    if variable not in scope:
        return factor
    idx = scope.index(variable)

    # new scope without variable.  
    new_scope = []
    for var in scope:
        if var != variable:
            new_scope.append(var)

    new_name = factor.name + "_sum_" + variable.name
    new_factor = Factor(new_name, new_scope)

    # new_assignments = []

    domain = list(map(Variable.domain, new_scope))

    assignments = get_assignments(domain, [], 0)

    new_assigns = []
    for assignment in assignments:
        count = 0
        for val in variable.domain():
            temp = list(assignment)
            temp.insert(idx, val)
            # print(temp)
            count += factor.get_value(temp)

        sum_out_assignment = list(assignment)
        sum_out_assignment.append(count)
        new_assigns.append(sum_out_assignment)

    new_factor.add_values(new_assigns)
    return new_factor

    
    


def multiply(factor_list):
    '''
    Multiply a list of factors together.
    Do not modify any of the input factors. 

    :param factor_list: a list of Factor objects.
    :return: a new Factor object resulting from multiplying all the factors in factor_list.
    '''
    if len(factor_list) == 1:
        return factor_list[0]
    

    new_factor = factor_list[0]
    for i in range(len(factor_list) - 1):
        new_factor = multiply_two_factors(new_factor, factor_list[i+1])
    
    return new_factor


def multiply_two_factors(factor1, factor2):

    new_scope = []
    new_values =[]
    
    for var in factor1.get_scope():
        if var not in new_scope:
            new_scope.append(var)

    for var in factor2.get_scope():
        if var not in new_scope:
            new_scope.append(var)


    new_name = "Mult[" + factor1.name + ", " + factor2.name + "]"
    new_factor = Factor(new_name, new_scope)

    domain = list(map(Variable.domain, new_scope))
    assignments = get_assignments(domain, [])

    for assignment in assignments:
        idx1 = []
        idx2 = []
        for var in factor1.get_scope():
            idx1.append(new_scope.index(var))
        for var in factor2.get_scope():
            idx2.append(new_scope.index(var))

        val1 = factor1.get_value([assignment[i] for i in idx1])
        val2 = factor2.get_value([assignment[i] for i in idx2])

        temp = list(assignment)
        temp.append(val1*val2)

        new_values.append(temp)

    new_factor.add_values(new_values)
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
        Pr(A='b'|B=1, C='c') = 0.24, and 
        Pr(A='c'|B=1, C='c') = 0.26.

    '''
    restricted_factors = []
    const = []

    for factor in bayes_net.factors():
        new_factor = factor
        
        for var in EvidenceVars:
            if var in new_factor.get_scope():
                new_factor = restrict(new_factor, var, var.get_evidence())
        
        if new_factor.get_scope() != []:
            restricted_factors.append(new_factor)
        else:
            const.append(new_factor)

    # for factor in restricted_factors:
    #     print(factor.name)
    #     factor.print_table()

    hidden = []
    for variable in bayes_net.variables():
        if variable not in EvidenceVars and variable != var_query:
            hidden.append(variable)

    remaining_factors = []

    for var in hidden:
        # print(var.name)
        mul_factors = []
        for factor in restricted_factors:
            if var in factor.get_scope():
                mul_factors.append(factor)

    
        if mul_factors == []:
            continue

        new_factor = multiply(mul_factors)


        new_factor = sum_out(new_factor, var)

        restricted_factors.append(new_factor)
        

        for factor in mul_factors:
            restricted_factors.remove(factor)





    # for factor in restricted_factors:
    #     factor.print_table()

    for factor in restricted_factors:
        if var_query in factor.get_scope():
            remaining_factors.append(factor)
        

    final_factor = multiply(remaining_factors)

        
        
    final_factor = normalize(final_factor)


    return final_factor




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
    ### READ IN THE DATA
    # counter = 0
    input_data = []
    with open(data_file, newline='') as csvfile:
        reader = csv.reader(csvfile)
        headers = next(reader, None) #skip header row
        for row in reader:
            # counter += 1
            input_data.append(row)

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
    
    #class_var = Variable("Salary", ['<50K', '>=50K'])
    var_list = {}
    for key in variable_domains:
        var_list[key] = Variable(key, variable_domains[key])

    factors = []

    class_values = {class_var.domain()[0]: 0, class_var.domain()[1]: 0}

    feature_value = {
        class_var.domain()[0]: {feature: {} for feature in variable_domains if feature != "Salary"},
        class_var.domain()[1]: {feature: {} for feature in variable_domains if feature != "Salary"}
    }
    # print(feature_value)
    # print(feature_value["<50K"]["Work"])
    count = 0
    for row in input_data:
        # print(row)
        class_values[row[-1]] += 1
        for key in variable_domains:
            if key == "Salary":
                continue

            if row[headers.index(key)] not in feature_value[row[-1]][key]:
                feature_value[row[-1]][key][row[headers.index(key)]] = 1
            else:
                feature_value[row[-1]][key][row[headers.index(key)]] += 1
        # for i in range(len(row) - 1):
        #     if row[i] == "Professional":
        #         count += 1
        #     if row[i] not in feature_value[row[-1]]:
        #         feature_value[row[-1]][row[i]] = 1
        #     else:
        #         feature_value[row[-1]][row[i]] += 1




    total_individuals = len(input_data) - 1
    # print(total_individuals)
    class_probability = {k: v / total_individuals for k, v in class_values.items()}
    # for k, v in class_values.items():
    #     print(v)
    #     print(total_individuals)
    #     print(v / total_individuals)


    class_factor = Factor(f"P({class_var.name})", [var_list[class_var.name]])
    f = []
    for val in class_var.domain():
        # print(class_probability[val])
        f.append([val, class_probability[val]])

    class_factor.add_values(f)

    # class_factor.print_table()
    factors.append(class_factor)

    
    for key in variable_domains:
        if key == class_var.name:
            continue
        feature_factor_values = []
        for feature_val in variable_domains[key]:
            for class_val in class_var.domain():
                if feature_val not in feature_value[class_val][key]:
                    prob = 0
                else:
                    #TODO: this is not correct!
                    #TODO: There is an error in counting, since there is a little overlap between variables. 
                    # Education has some Professionals, and occupation also has some professionals for example. 


                    # if key == "Education":
                    #     if feature_val == "Professional":
                    #         print(feature_value[class_val][feature_val])
                    #         print(feature_value[class_val][feature_val] / class_values[class_val])
                    prob = feature_value[class_val][key][feature_val] / class_values[class_val]
                    feature_factor_values.append([feature_val, class_val, prob])

                    
        feature_factor = Factor(f"P({key}|{class_var.name})", [var_list[key], var_list[class_var.name]])
        feature_factor.add_values(feature_factor_values)

        # feature_factor.print_table()

        factors.append(feature_factor)


    new_var_list = [var_list[key] for key in var_list]
    bayes_net = BN("NaiveBayesModel", new_var_list, factors)
    return bayes_net
        



def compute_probability(bayes_net, e1):
    # gender_var = bayes_net.get_variable("Gender")
    salary_var = bayes_net.get_variable("Salary")

    for evidence in e1:
        evidence[0].set_evidence(evidence[1])

    factor = ve(bayes_net, salary_var, [ev[0] for ev in e1])
    return factor.get_value([">=50K"])
    
    


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
    input_data = []
    with open('data/adult-test.csv', newline='') as csvfile:
        reader = csv.reader(csvfile)
        headers = next(reader, None) #skip header row
        for row in reader:
            # counter += 1
            input_data.append(row)

        
    if question not in [1, 2, 3, 4, 5, 6]:
        raise ValueError("Invalid question number.")
    
    if question == 1:
        e1_count = 0
        e2_count = 0
        for data in input_data:
            if data[headers.index("Gender")] == "Female":
                e1 = [(bayes_net.get_variable("Work"), data[headers.index("Work")]), 
                    (bayes_net.get_variable("Education"), data[headers.index("Education")]), 
                    (bayes_net.get_variable("Occupation"), data[headers.index("Occupation")]), 
                    (bayes_net.get_variable("Relationship"), data[headers.index("Relationship")])]
                
                e2 = [(bayes_net.get_variable("Work"), data[headers.index("Work")]), 
                    (bayes_net.get_variable("Education"), data[headers.index("Education")]), 
                    (bayes_net.get_variable("Occupation"), data[headers.index("Occupation")]), 
                    (bayes_net.get_variable("Relationship"), data[headers.index("Relationship")]),
                    (bayes_net.get_variable("Gender"), data[headers.index("Gender")])]
                
                prob1 = compute_probability(bayes_net, e1)
                prob2 = compute_probability(bayes_net, e2)
                if prob1 > prob2:
                    e1_count += 1
                else:
                    e2_count += 1
        

        e2_count = e1_count - e2_count

        return ((e2_count / e1_count) ) * 100
                  
    
    elif question == 2:
        e1_count = 0
        e2_count = 0
        for data in input_data:
            if data[headers.index("Gender")] == "Male":
                e1 = [(bayes_net.get_variable("Work"), data[headers.index("Work")]), 
                    (bayes_net.get_variable("Education"), data[headers.index("Education")]), 
                    (bayes_net.get_variable("Occupation"), data[headers.index("Occupation")]), 
                    (bayes_net.get_variable("Relationship"), data[headers.index("Relationship")])]
                
                e2 = [(bayes_net.get_variable("Work"), data[headers.index("Work")]), 
                    (bayes_net.get_variable("Education"), data[headers.index("Education")]), 
                    (bayes_net.get_variable("Occupation"), data[headers.index("Occupation")]), 
                    (bayes_net.get_variable("Relationship"), data[headers.index("Relationship")]),
                    (bayes_net.get_variable("Gender"), data[headers.index("Gender")])]
                
                prob1 = compute_probability(bayes_net, e1)
                prob2 = compute_probability(bayes_net, e2)
                if prob1 > prob2:
                    e1_count += 1
                else:
                    e2_count += 1
        
        return (e1_count / e2_count) * 100
    
    elif question == 3:
        actual_count = 0
        predicted_count = 0
        for data in input_data:
            if data[headers.index("Gender")] == "Female":
                e1 = [(bayes_net.get_variable("Work"), data[headers.index("Work")]), 
                    (bayes_net.get_variable("Education"), data[headers.index("Education")]), 
                    (bayes_net.get_variable("Occupation"), data[headers.index("Occupation")]), 
                    (bayes_net.get_variable("Relationship"), data[headers.index("Relationship")])]
                
                prob = compute_probability(bayes_net, e1)
                if prob > 0.5:
                    predicted_count += 1
                    if data[headers.index("Salary")] == ">=50K":
                        actual_count += 1
        # print(predicted_count)
        # print(actual_count)
        return (actual_count / predicted_count) * 100
                

  
            

    elif question == 4:
        actual_count = 0
        predicted_count = 0
        for data in input_data:
            if data[headers.index("Gender")] == "Male":
                e1 = [(bayes_net.get_variable("Work"), data[headers.index("Work")]), 
                    (bayes_net.get_variable("Education"), data[headers.index("Education")]), 
                    (bayes_net.get_variable("Occupation"), data[headers.index("Occupation")]), 
                    (bayes_net.get_variable("Relationship"), data[headers.index("Relationship")])]
                
                prob = compute_probability(bayes_net, e1)
                if prob > 0.5:
                    predicted_count += 1
                    if data[headers.index("Salary")] == ">=50K":
                        actual_count += 1

        return (actual_count / predicted_count) * 100

    elif question == 5:
        actual_count = 0
        predicted_count = 0
        for data in input_data:
            if data[headers.index("Gender")] == "Female":
                e1 = [(bayes_net.get_variable("Work"), data[headers.index("Work")]), 
                    (bayes_net.get_variable("Education"), data[headers.index("Education")]), 
                    (bayes_net.get_variable("Occupation"), data[headers.index("Occupation")]), 
                    (bayes_net.get_variable("Relationship"), data[headers.index("Relationship")]), ]
                
                prob = compute_probability(bayes_net, e1)
                if prob > 0.5:
                    predicted_count += 1

                actual_count += 1

        return (predicted_count / actual_count) * 100

    elif question == 6:
        actual_count = 0
        predicted_count = 0
        for data in input_data:
            if data[headers.index("Gender")] == "Male":
                e1 = [(bayes_net.get_variable("Work"), data[headers.index("Work")]), 
                    (bayes_net.get_variable("Education"), data[headers.index("Education")]), 
                    (bayes_net.get_variable("Occupation"), data[headers.index("Occupation")]), 
                    (bayes_net.get_variable("Relationship"), data[headers.index("Relationship")])]
                
                prob = compute_probability(bayes_net, e1)
                if prob > 0.5:
                    predicted_count += 1
                actual_count += 1

        return (predicted_count / actual_count) * 100


if __name__ == '__main__':
    nb = naive_bayes_model('data/adult-train.csv')
    temp = []

    for variable in nb.variables():
        if variable.name == "Salary":
            x = variable

        # if variable.name == "Gender":
        #     temp.append(variable)
        #     variable.set_evidence("Female")
        
    #     if variable.name == "Work":
    #         temp.append(variable)
    #         variable.set_evidence("Private")

    #     if variable.name == "Education":
    #         temp.append(variable)
    #         variable.set_evidence("HS-Graduate")

    #     if variable.name == "MaritalStatus":
    #         temp.append(variable)
    #         variable.set_evidence("Married")

    #     if variable.name == "Occupation":
    #         temp.append(variable)
    #         variable.set_evidence("Manual Labour")

    #     if variable.name == "Relationship":
    #         temp.append(variable)
    #         variable.set_evidence("Husband")

    #     if variable.name == "Race":
    #         temp.append(variable)
    #         variable.set_evidence("White")

    #     if variable.name == "Country":
    #         temp.append(variable)
    #         variable.set_evidence("North-America")

        

    

    print("here")    

    # ve(nb, x, []).print_table()
    # print(22654 / 30160)



    for i in range(1,7):
        print("explore(nb,{}) = {}".format(i, explore(nb, i)))
