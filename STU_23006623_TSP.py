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
    n = len(gene_pool)
    
    if n == 0:
        return gene_pool
    # -1 because the random.randint method is inclusive
    elite_range = n - 1 if elite_threshold == 0 else int((n - 1) * elite_threshold)
    
    elites = gene_pool[:elite_range + 1]
    
    progeny = [best_solution]
    
    while len(progeny) < n * 3:
        parent1 = random.choice(elites)
        parent2 = random.choice(gene_pool)
        
        bred_offspring = breed(parent1, parent2)
        mutated_offspring = mutate(bred_offspring, mutation_rate)
        
        progeny.append(mutated_offspring)
        progeny.append(bred_offspring)

    sorted_progeny = sorted(progeny, key=lambda x: calculate_path(x))
    
    return sorted_progeny[:n]


def breed(parent1: list, parent2: list) -> list:
    offspring_p1 = []
    offspring_p2 = []
    
    gene_a = int(random.random() * len(parent1))
    gene_b = int(random.random() * len(parent1))
    
    start_pt = min(gene_a, gene_b)
    end_pt = max(gene_a, gene_b)

    for i in range(start_pt, end_pt):
        offspring_p1.append(parent1[i])
        
    offspring_p2 = [item for item in parent2 if item not in offspring_p1]

    offspring = offspring_p1 + offspring_p2
    return offspring


def mutate(child: list, mutation_rate: float) -> list:    
    for _ in range(len(child)):
        if random.random() < mutation_rate:
            idx1 = int(random.random() * len(child))
            idx2 = int(random.random() * len(child))

            child[idx1], child[idx2] = child[idx2], child[idx1]

    return child


def iterator(gene_pool: list, iterations: int, mutation_rate: float, elite_threshold: float) -> tuple:
    best_solution = None
    current_gene_pool = gene_pool
    
    for _ in range(iterations):
        current_gene_pool, best_solution = fitness_function(current_gene_pool, best_solution)
        current_gene_pool = mating_function(current_gene_pool, best_solution, mutation_rate, elite_threshold)
    
    return best_solution