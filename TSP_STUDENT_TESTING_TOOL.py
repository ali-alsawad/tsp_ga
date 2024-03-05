import random
import math
import csv
import turtle
import copy
import inspect
import timeit

from inspect import signature

provided_map = []
    
def testing_procedure (current_student):
    global student
    global current_test
    
    student = current_student

    setup_feedback_array(12)

    # Test 1 - Student details function implementation
    start = timeit.default_timer()
    current_test = 0
    test_student_details("Test 1 - Student details function implementation")
    stop = timeit.default_timer()
    feedback_generator("")

    # Test 2 - Map generation function implementation
    start = timeit.default_timer()
    current_test = 1
    test_map_generation(300, 300, 10, "Test 2 - Map generation function implementation")
    stop = timeit.default_timer()
    feedback_generator("")

    # Test 3 - Path calculation function implementation
    start = timeit.default_timer()
    current_test = 2
    test_path_calculation(provided_map, "Test 3 - Path calculation function implementation")
    stop = timeit.default_timer()
    feedback_generator("")

    # Test 4 - Print function implementation
    start = timeit.default_timer()
    current_test = 3
    test_print_function(0, "red", 4, "green" ,1 ,provided_map, "Test 4 - Print function implementation")
    stop = timeit.default_timer()
    feedback_generator("")

    # Test 5 - Nearest neighbour algorithm function implementation
    start = timeit.default_timer()
    current_test = 4
    test_algorithm(provided_map,"Test 5 - Nearest neighbour algorithm function implementation")
    stop = timeit.default_timer()
    feedback_generator("")
    
    # Test 6 - Test gene pool creation
    start = timeit.default_timer()
    current_test = 5
    test_create_population(49, provided_map, "Test 6 - Create population function implementation")
    stop = timeit.default_timer()
    feedback_generator("")

    # Test 7 - Test fitness function
    start = timeit.default_timer()
    current_test = 6
    test_fitness_function(49, provided_map, "Test 7 - Fitness function implementation")
    stop = timeit.default_timer()
    feedback_generator("")
    
    # Test 8 - Test mating function
    start = timeit.default_timer()
    current_test = 7
    test_mating_function(49, provided_map, 0.5, 0.1, "Test 8 - Mating function implementation")
    stop = timeit.default_timer()
    feedback_generator("")

    # Test 9 - Test breeding function
    start = timeit.default_timer()
    current_test = 8
    test_breed(49, provided_map, "Test 9 - Breeding function implementation")
    stop = timeit.default_timer()
    feedback_generator("")
    
    # Test 10 - Test mutate function
    start = timeit.default_timer()
    current_test = 9
    test_mutate(provided_map, 0.5, "Test 10 - Mutation function implementation")
    stop = timeit.default_timer()
    feedback_generator("")
    

    # Test 11 - Test iterator function
    start = timeit.default_timer()
    current_test = 10
    test_iterator(49, 12, 0.01, 0.1, provided_map, "Test 11 - Iterator function implementation")
    stop = timeit.default_timer()
    feedback_generator("")

    # Test 12 - Test full genetic algorithm
    start = timeit.default_timer()
    current_test = 11
    population = 5
    iterations = 20
    mutation_rate = 0.01
    elite_threshold = 0.1
    test_genetic_algorithm(provided_map, population, iterations, mutation_rate, elite_threshold, f"Test 12 - Full run of genetic algorithm with population of {population}, {iterations} iterations, a mutation rate of {mutation_rate} and an elite threshold of {elite_threshold}")
    stop = timeit.default_timer()
    feedback_generator("")
    print("")

    # final report
    for x in feedback_array:
        for i in x:
            print(i[0])

########### grading tools ################

def setup_feedback_array(number_of_tests):
    global feedback_array
    feedback_array = []
    for x in range(number_of_tests):
        feedback_array.append([])

def feedback_generator(feedback):
    feedback_array[current_test].append([f"{feedback}"])

 
###################################################################################

def student_testing(test_item,test_item_str,arguments):
    function_exists = (test_item_str in dir(student))
    if function_exists == True:
        feedback_generator("Function found successfully")
        if len(signature(test_item).parameters) == arguments:
            feedback_generator("Correct number of arguments found")
        else:
            feedback_generator("Incorrect number of arguments found")

########################### student details testing function ####################################
def test_student_details(testID):
    global student
    feedback_generator(testID)
    print(testID)

    student_id = None
    student_username = None
    
    # Test inputs
    try:
        student_testing(student.student_details,"student_details",0)
    except:
        feedback_generator("student_details function not found")
        
    # Test outputs
    try:
        student_id, student_username = student.student_details()
        print("student_id =",student_id)
        print("student_username =",student_username)
        feedback_generator("Correct number of values returned")
    except:
        feedback_generator("Incorrect number of values returned, should have been 2 values returned")

    # check student ID returned from function
    if student_id is not None:
        if type(student_id) == int:
            feedback_generator("Student ID is using the correct data type (int)")
        else:
            feedback_generator("Student ID is using the incorrect data type, it should be (int)")
    else:
        feedback_generator("Cannot access student ID details")

    # check student username returned from function
    if student_username is not None:
        if type(student_username) == str:
            feedback_generator("Student username is using the correct data type (str)")
        else:
            feedback_generator("Student username is using the incorrect data type, should be (str)")
    else:
        feedback_generator("Cannot access student username details")


########################### map generation testing function ##################################
def test_map_generation(xrange,yrange,cities,testID):
    feedback_generator(testID)
    print(testID)
    
    generated_map = None

    # Test inputs
    try:
        student_testing(student.generate_map,"generate_map",3)
    except:
        feedback_generator("generate_map function not found")

    # Test outputs
    try:
        generated_map = student.generate_map(xrange,yrange,cities)
        print("generated_map =",generated_map)
        feedback_generator("Correct number of values returned")
    except:
        feedback_generator("Incorrect number of values returned, should have been 1 value returned")

    # Test map is not empty
    if generated_map is not None:

        # Test map is correct data type
        if type(generated_map) == list:
            feedback_generator("Correct data type for map generated")

            # Test correct number of cities
            if len(generated_map) == cities:
                feedback_generator("Map generated contained correct number of cities")
                correct_ranges = True
                positive_range = 0
                negative_range = 0
                for x in generated_map:
                    if x[0] > xrange or x[0] < -xrange:
                        correct_ranges = False
                    if x[1] > yrange or x[1] < -yrange:
                        correct_ranges = False

                    if correct_ranges == True:
                        if x[0] < 0:
                            negative_range += 1
                        else:
                            positive_range += 1
                        if x[1] < 0:
                            negative_range += 1
                        else:
                            positive_range += 1

                # Test cities within ranges
                if correct_ranges == True:
                    feedback_generator("Map generated coordinates within defined parameters")
                else:
                    feedback_generator("Map generated coordinates outside of defined parameters")

                # Test cities using full range
                if (positive_range > 0) and (negative_range >0):
                    feedback_generator("Map generated coordinates utilising full positive and negative range")
                else:
                    if(positive_range == 0):
                        feedback_generator("Map generated coordinates not utilising positive range")
                    else:
                        feedback_generator("Map generated coordinates not utilising negative range")
            else:
                feedback_generator("Map generated contained incorrect number of cities")
        else:
            feedback_generator("Incorrect data type for map generated")
    else:
        feedback_generator("Cannot access generated map")

########################### path calculation testing function ##################################
def test_path_calculation(provided_map,testID):
    feedback_generator(testID)
    print(testID)
    
    selected_map = copy.deepcopy(provided_map)
    selected_map_path_calculation = None

    # Test inputs
    try:
        student_testing(student.calculate_path,"calculate_path",1)
    except:
        feedback_generator("calculate_path function not found")

    # Test outputs
    try:
        selected_map_path_calculation = round(student.calculate_path(selected_map),4)
        print("selected_map_path_calculation =",selected_map_path_calculation)
        feedback_generator("Correct number of values returned")
    except:
        feedback_generator("Incorrect number of values returned, should have been 1 value returned")

    # Test path calculation is not empty
    if selected_map_path_calculation is not None:

        # Test data type of path calculation output
        if (type(selected_map_path_calculation) == float):
            feedback_generator("Correct data type used to output path calculation")
        else:
            feedback_generator("Incorrect data type used to output path calculation")
    else:
        feedback_generator("Cannot access path calculation")


########################### print testing function ##################################
def test_print_function(speed,s_color,s_thickness,m_color,m_thickness,provided_map,testID):
    feedback_generator(testID)
    print(testID)
    
    selected_map = copy.deepcopy(provided_map)

    # Test inputs
    try:
        student_testing(student.print_map,"print_map",4)
    except:
        feedback_generator("print_map function not found")
        
    # Test outputs
    try:
        student.print_map(speed,s_color,s_thickness,selected_map)
        feedback_generator("printing map function accessed")
    except:
        feedback_generator("Cannot access print map function")

    # MANUAL CHECK REQUIRED!!! 

########################### algorithm testing function ##################################
def test_algorithm(provided_map,testID):
    feedback_generator(testID)
    print(testID)
    
    selected_map = copy.deepcopy(provided_map)
    nearest_neighbour_map_optimised = None

    # Test inputs
    try:
        student_testing(student.nearest_neighbour_algorithm, "nearest_neighbour_algorithm", 1)  
    except:
        feedback_generator("Nearest neighbour algorithm not found")

    # Test outputs
    try:
        nearest_neighbour_map_optimised = student.nearest_neighbour_algorithm(selected_map)
        print("nearest_neighbour_map_optimised =",nearest_neighbour_map_optimised)
        feedback_generator("Correct number of values returned")
    except:
        feedback_generator("Incorrect number of values returned, should have been 1 value returned")

    # Test map is not empty
    if nearest_neighbour_map_optimised is not None:

        # Test map is correct data type
        if type(nearest_neighbour_map_optimised) == list:
            feedback_generator("Returning correct data type (list)")

            # Test map is using correct locations
            check_locations = 0   
            for x in provided_map:
                for i in nearest_neighbour_map_optimised:
                    if x == i:
                        check_locations +=1
            if len(provided_map) == check_locations:
                feedback_generator(f"Nearest neighbour algorithm is correct, using all original provided locations")
            else:
                feedback_generator("Nearest neighbour algorithm is incorrect and is not using all original provided locations")
        else:
            feedback_generator("Returning incorrect data type, should have been list")
    else:
        feedback_generator("Cannot access nearest neighbour algorithm")

##################### testing population creation ########################################
def test_create_population(population, provided_map, testID):
    feedback_generator(testID)
    print(testID)
    
    selected_map = copy.deepcopy(provided_map)
    gene_pool = None

    # Test inputs
    try:
        student_testing(student.create_population, "create_population", 2)  
    except:
        feedback_generator("Genetic algorithm - create_population function not found")
    
    # Test outputs
    try:
        gene_pool = student.create_population(population, selected_map)
        print("gene_pool = ",gene_pool)
        feedback_generator("Correct number of values returned")
    except:
        feedback_generator("Incorrect number of values returned, should have been 1 value returned")

    # Test gene pool is not empty
    if gene_pool is not None:

        # Test size of gene pool
        if len(gene_pool) == population:
            feedback_generator("Gene pool contains correct number of individuals")
        else:
            feedback_generator("Gene pool contains incorrect number of individuals")

        # Test validity of gene pool
        if marker_check_gene_pool_validity(gene_pool, provided_map):
            feedback_generator("Individuals in gene pool are correct and contain all locations")
        else:
            feedback_generator("Individuals in gene pool contain location errors (not all destinations match the map)")
        
        # Test randomness of gene pool
        shuffle_error = False
        for x in gene_pool:
            if x == provided_map:
                shuffle_error = True
        if shuffle_error == False:
            feedback_generator("Individuals in gene pool are in a random order")
        else:
            feedback_generator("Individuals in gene pool are not in a random order, this is incorrect")
    else:
        feedback_generator("Cannot access Genetic algorithm - create_population function")

####################### testing fitness function #################################
def test_fitness_function(population, provided_map, testID):
    feedback_generator(testID)
    print(testID)

    provided_gene_pool = student_create_population(population, provided_map)

    gene_pool_selected = copy.deepcopy(provided_gene_pool)
    best_solution_selected = provided_map
    provided_best_solution = provided_map

    gene_pool = None
    best_solution = None

    # Test inputs
    try:
        student_testing(student.fitness_function, "fitness_function", 2)  
    except:
        feedback_generator("Genetic algorithm - fitness_function not found")
    
    # Test outputs
    try:
        gene_pool, best_solution = student.fitness_function(gene_pool_selected, best_solution_selected)
        print("gene_pool =",gene_pool)
        print("best_solution =",best_solution)
        feedback_generator("Correct number of values returned")
    except:
        feedback_generator("Incorrect number of values returned, should have been 2 values returned")

    # Test gene pool is not empty
    if gene_pool is not None:

        # Test size of gene pool
        if len(gene_pool) == len(provided_gene_pool):
            feedback_generator("Gene pool contains correct number of individuals")
        else:
            feedback_generator("Gene pool contains incorrect number of individuals")

        # Test validity of gene pool
        if marker_check_gene_pool_validity(gene_pool, provided_map):
            feedback_generator("Individuals in gene pool are correct and contain all locations")
        else:
            feedback_generator("Individuals in gene pool contain location errors (not all destinations match the map)")

    else:
        feedback_generator("Cannot access Genetic algorithm - fitness_function")

######################## testing mating function ######################################
def test_mating_function(population, provided_map, mutation_rate, elite_threshold, testID):
    feedback_generator(testID)
    print(testID)

    provided_gene_pool = student_create_population(population, provided_map)
    gene_pool = copy.deepcopy(provided_gene_pool)
    best_solution = provided_map
    provided_best_solution = provided_map

    gene_pool_1 = None
    gene_pool_2 = None

    # Test inputs
    try:
        student_testing(student.mating_function, "mating_function", 4)
    except:
        feedback_generator("Genetic algorithm - mating_function not found")
    
    # Test outputs
    try:
        gene_pool_1 = student.mating_function(gene_pool, best_solution, mutation_rate, elite_threshold)
        gene_pool_2 = student.mating_function(gene_pool, best_solution, mutation_rate, elite_threshold)
        print("gene_pool_1 =",gene_pool_1)
        print("gene_pool_2 =",gene_pool_2)
        feedback_generator("Correct number of values returned")
    except:
        feedback_generator("Incorrect number of values returned, should have been 1 value returned")

    # Test gene pool is not empty
    if gene_pool_1 is not None and gene_pool_2 is not None:

        # Test gene pool contains correct population
        if len(gene_pool_1) == len(provided_gene_pool):
            feedback_generator("Gene pool contains correct number of individuals")
        else:
            feedback_generator("Gene pool contains incorrect number of individuals")

        # Test gene pool validity
        if marker_check_gene_pool_validity(gene_pool_1, provided_map):
            feedback_generator("Individuals in gene pool are correct and contain all locations")
        else:
            feedback_generator("Individuals in gene pool contain location errors (not all destinations match the map)")

        # Test mutations have taken place
        if gene_pool_1 != gene_pool_2 and gene_pool != gene_pool_1:
            feedback_generator("Your mating function appears to be working")
        else:
            feedback_generator("Your mating function appears to be functioning incorrectly")
    else:
        feedback_generator("Cannot access Genetic algorithm - mating_function")

######################### testing breed function #####################################
def test_breed(population, provided_map, testID):
    feedback_generator(testID)
    print(testID)
    provided_gene_pool = student_create_population(population, provided_map)

    parent_1 = provided_map
    parent_2 = provided_map
    parent_1_copy = copy.deepcopy(parent_1)
    parent_2_copy = copy.deepcopy(parent_2)

    child = None

    # Test inputs
    try:
        student_testing(student.breed, "breed", 2)
    except:
        feedback_generator("Genetic algorithm - breed function not found")
    
    # Test outputs
    try:
        child = student.breed(parent_1, parent_2)
        print("child =",child)
        feedback_generator("Correct number of values returned")
    except:
        feedback_generator("Incorrect number of values returned, should have been 1 value returned")

    # Test child is not empty
    if child is not None:

        # Test child is correct size
        if len(child) == len(provided_map):
            feedback_generator("Individual output contains correct number of locations")
        else:
            feedback_generator("Individual output contains incorrect number of locations")

        # Test all locations are present
        locations = 0
        for x in child:
            locations += provided_map.count(x)
        if locations == len(provided_map):
            feedback_generator("Individual output is correct and contain all locations")

            # Test genetic crossover is working correctly
            checking_dna = False
            attempts = 10
            for x in range(attempts):
                if checking_dna == False:
                    child = student.breed(parent_1, parent_2)
                    checking_dna = marker_check_parents(child, parent_1_copy, parent_2_copy, provided_map)
            if checking_dna == True:
                feedback_generator("Breed function is working correctly")
            else:
                feedback_generator("Breed function is not working correctly")
        else:
            feedback_generator("individual output contains location errors (not all destinations match the map)")
    else:
        feedback_generator("Cannot access Genetic algorithm - breed function")

########################### testing mutate function ######################################
def test_mutate(provided_map, mutation_rate, testID):
    feedback_generator(testID)
    print(testID)
    selected_map = copy.deepcopy(provided_map)

    child = None

    # Test inputs
    try:
        student_testing(student.mutate, "mutate", 2)
    except:
        feedback_generator("Genetic algorithm - mutate function not found")
    
    # Test outputs
    try:
        child = student.mutate(selected_map, mutation_rate)
        print("child =",child)
        feedback_generator("Correct number of values returned")
    except:
        feedback_generator("Incorrect number of values returned, should have been 1 value returned")

    # Test child is not empty
    if child is not None:

        # Test child is correct size
        if len(child) == len(provided_map):
            feedback_generator("Individual output contains correct number of locations")
        else:
            feedback_generator("Individual output contains incorrect number of locations")

        # Test all locations are within child
        locations = 0
        for x in child:
            locations += provided_map.count(x)
        if locations == len(provided_map):
            feedback_generator("Individual output is correct and contain all locations")
        else:
            feedback_generator("Individual output contains location errors (not all destinations match the map)")

        # Test for mutations taking place
        if child != provided_map:
            feedback_generator("Individual output has changed suggesting that some mutations have taken place")
            count_mutations = 0
            for x in range(len(child)-1):
                if child[x] != provided_map[x]:
                    count_mutations += 1
            feedback_generator(f"{count_mutations} mutations have taken place with a mutation rate of {mutation_rate}")
        else:
            feedback_generator("Individual output has not changed suggesting that no mutations have taken place and there is an error with the function")
    else:
        feedback_generator("Cannot access Genetic algorithm - mutate function")
        
########################### testing iterator function ######################################
def test_iterator(population, iterations, mutation_rate, elite_threshold, provided_map, testID):
    feedback_generator(testID)
    print(testID)
    provided_gene_pool = student_create_population(population, provided_map)

    best_solution = None

    # Test inputs
    try:
        student_testing(student.iterator, "iterator", 4)
    except:
        feedback_generator("Genetic algorithm - iterator function not found")
    
    # Test outputs
    try:
        best_solution = student.iterator(provided_gene_pool, iterations, mutation_rate, elite_threshold)
        print("best_solution =",best_solution)
        feedback_generator("Correct number of values returned")
    except:
        feedback_generator("Incorrect number of values returned, should have been 1 value returned")

    # Test best_solution is not empty
    if best_solution is not None:    

        # Test best_solution is correct size
        if len(best_solution) == len(provided_map):
            feedback_generator("Individual output contains correct number of locations")
        else:
            feedback_generator("Individual output contains incorrect number of locations")

        # Test all locations are within best solution
        locations = 0
        for x in best_solution:
            locations += provided_map.count(x)
        if locations == len(provided_map):
            feedback_generator("Individual output is correct and contain all locations")
        else:
            feedback_generator("Individual output contains location errors (not all destinations match the map)")
    else:
        feedback_generator("Cannot access Genetic algorithm - iterator function")
        
############################ testing full genetic algorithm function #############################
def test_genetic_algorithm(provided_map, population, iterations, mutation_rate, elite_threshold, testID):
    feedback_generator(testID)
    print(testID)
    selected_map = copy.deepcopy(provided_map)
    best_solution = None

    # Test inputs
    try:
        student_testing(student.genetic_algorithm, "genetic_algorithm", 5)
    except:
        feedback_generator("Genetic algorithm not found")
    
    # Test outputs
    try:
        best_solution = student.genetic_algorithm(selected_map, population, iterations, mutation_rate, elite_threshold)
        print("best_solution =",best_solution)
        feedback_generator("Correct number of values returned")
    except:
        feedback_generator("Incorrect number of values returned, should have been 1 value returned")

    # Test map is not empty
    if best_solution is not None:

        # Test map is correct data type
        if type(best_solution) == list:
            feedback_generator("Returning correct data type (list)")

            # Test map is using correct locations
            check_locations = 0   
            for x in provided_map:
                for i in best_solution:
                    if x == i:
                        check_locations +=1
            if len(provided_map) == check_locations:
                feedback_generator(f"Genetic algorithm is correct, using all original provided locations ")
            else:
                feedback_generator("Genetic algorithm is incorrect and is not using all original provided locations")
        else:
            feedback_generator("Returning incorrect data type, should have been list")
    else:
        feedback_generator("Cannot access genetic algorithm")

#############################################################################
def marker_check_parents(child, parent_1, parent_2, provided_map):
    parent_1_genes = 0
    parent_2_genes = 0
    for x in range(len(child)):
        if child[x] == parent_1[x]:
            parent_1_genes += 1
        if child[x] == parent_2[x]:
            parent_2_genes += 1
    parents_dna_used = False
    if parent_1_genes > 0 or parent_2_genes > 0:
        parents_dna_used = True
    if child == parent_1 or child == parent_2 or child == provided_map:
        parents_dna_used = False
    return parents_dna_used
        
    
    
def marker_check_gene_pool_validity(gene_pool, provided_map):
    validated_map = True
    location_errors = 0
    for x in gene_pool:
        temp_locations = 0
        for i in provided_map:
            temp_locations += x.count(i)
        if temp_locations != len(provided_map):
            location_errors += 1
    if location_errors > 0:
        validated_map = False

    return validated_map

def student_create_population(population, provided_map):
    
    gene_pool = None
    try:
        gene_pool = student.create_population(population, provided_map)
        feedback_generator("Accessing student create_population funtion to test this function")
    except:
        feedback_generator("Unable to access student create_population funtion so cannot perform this test")

    return gene_pool 

###########################################################################################################################

def open_map():
    with open('TSP.csv', newline='', encoding='utf-8') as f:
        reader = csv.reader(f, delimiter=',')
        rows = [[int(row[0]), int(row[1])] for row in reader]
    for row in rows:
        provided_map.append(row)
    
open_map()


###########################################################################################################################


import STU_23006623_TSP
testing_procedure(STU_23006623_TSP)                          # this is calling the testing procedue with your file

