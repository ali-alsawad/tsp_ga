import random
import math
import csv
import turtle
import copy
import inspect
import timeit
#### code below this point do not add any libraries ####

def student_details() -> tuple:
    return 23006623, 'aa23aqe'


def generate_map(xrange: int, yrange: int, cities: int) -> list:
    return [[random.randint(-xrange, xrange), random.randint(-yrange, yrange)]
            for _ in range(cities)]


def calculate_distance(x1, y1, x2, y2) -> float:
    return math.sqrt((x2 - x1)**2 + (y2 - y1)**2)


def calculate_path(selected_map: list) -> float:
    calculated_path = 0.0

    for i in range(len(selected_map)):
        pos = i + 1
        if pos == len(selected_map):
            # if we are at the end of the list, we need to go back to the start
            pos = 0

        st_pt = selected_map[i]
        end_pt = selected_map[pos]

        calculated_path += calculate_distance(
            st_pt[0], st_pt[1], end_pt[0], end_pt[1])

    return calculated_path


def print_map(speed: int, color: str, thickness: int, selected_map: list) -> None:
    t = turtle.Turtle()
    t.speed(speed)
    t.color(color)
    t.pensize(thickness)

    start_pos = selected_map[0]
    t.penup()
    t.goto(start_pos)
    t.pendown()

    for i in range(1, len(selected_map)):
        t.goto(selected_map[i])
    t.goto(start_pos) 

    t.hideturtle()


def new_distances(st: list, points: list) -> list:
    return [calculate_distance(st[0], st[1], p[0], p[1]) for p in points]


def swap_last_pop(points: list, idx: int) -> list:
    points[idx], points[-1] = points[-1], points[idx]
    return points.pop()


def nearest_neighbour_algorithm(selected_map: list) -> list:
    if len(selected_map) == 0:
        return []
    
    temp = selected_map[:]
    # retaining order is not necessary so we can swap and pop for an O(1) operation
    optimized_map = [swap_last_pop(temp, 0)]


    while len(temp) > 0:
        # get distances from the latest point to all other remaining points
        updated_distances = new_distances(optimized_map[-1], temp)
        
        min_dist_idx = None
        for i, dist in enumerate(updated_distances):
            if min_dist_idx is None or dist < updated_distances[min_dist_idx]:
                min_dist_idx = i

        nearest_pt = swap_last_pop(temp, min_dist_idx)
        optimized_map.append(nearest_pt)

    return optimized_map


def genetic_algorithm(selected_map: list, population: int, iterations: int, mutation_rate: float, elite_threshold: float) -> list:
    gene_pool = create_population(population, selected_map)
    best_solution = iterator(gene_pool, iterations, mutation_rate, elite_threshold)
    return best_solution


def create_population(population: int, selected_map: list) -> list:
    gene_pool = []
    
    for _ in range(population):
        sample = random.sample(selected_map, len(selected_map))
        gene_pool.append(sample)
        
    return gene_pool


def fitness_function(gene_pool: list, best_solution: list) -> tuple:
    best_distance = float('inf')
    if best_solution is not None:
        best_distance = calculate_path(best_solution)

    distances = [calculate_path(gene) for gene in gene_pool]

    gene_distances = sorted(zip(gene_pool, distances), key=lambda x: x[1])

    sorted_gene_pool, sorted_distances = zip(*gene_distances)

    best_candidate = sorted_gene_pool[0]
    candidate_distance = sorted_distances[0]

    if candidate_distance < best_distance:
        best_solution = best_candidate

    return list(sorted_gene_pool), best_solution


def mating_function(gene_pool: list, best_solution: list, mutation_rate: float, elite_threshold: float) -> list:
    progeny = []
    sorted_gene_pool, _  = fitness_function(gene_pool, best_solution)
    n = len(sorted_gene_pool)
    # -1 because the random.randint method is inclusive
    elite_range = n - 1 if elite_threshold == 0 else int((n - 1) * elite_threshold)
    
    while len(progeny) < n:
        elite_parent_idx = random.randint(0, elite_range)
        normal_parent_idx = random.randint(0, n - 1)
        
        if n > 1:
            while elite_parent_idx == normal_parent_idx:
                normal_parent_idx = random.randint(0, n - 1)
            
        elite_parent = sorted_gene_pool[elite_parent_idx]
        normal_parent = sorted_gene_pool[normal_parent_idx]
        
        offspring = breed(normal_parent, elite_parent)
        offspring_gene = mutate(offspring, mutation_rate)

        progeny.append(offspring_gene)

    return progeny


def breed(parent1: list, parent2: list) -> list:
    n = len(parent1)
    offspring = []
    selected_genes = set()

    while len(selected_genes) < n:
        rand_pos = random.randint(0, n - 1)
        selected_gene = parent1[rand_pos] if random.random(
        ) < 0.5 else parent2[rand_pos]

        selected_gene_str = f'{selected_gene[0]}{selected_gene[1]}'
        if selected_gene_str not in selected_genes:
            selected_genes.add(selected_gene_str)
            offspring.append(selected_gene)

    return offspring


def mutate(child: list, mutation_rate: float) -> list:
    n = len(child)

    for i in range(n):
        if random.random() < mutation_rate:
            rand_pos = random.randint(0, n - 1)
            child[i], child[rand_pos] = child[rand_pos], child[i]

    return child


def iterator(gene_pool: list, iterations: int, mutation_rate: float, elite_threshold: float) -> tuple:
    best_solution = None
    current_gene_pool = gene_pool
    
    for _ in range(iterations):
        current_gene_pool, best_solution = fitness_function(current_gene_pool, best_solution)
        current_gene_pool = mating_function(current_gene_pool, best_solution, mutation_rate, elite_threshold)
    
    return best_solution