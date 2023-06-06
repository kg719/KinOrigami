#include <iostream>
#include <map>
#include <vector>
#include <cmath>
#include <queue>
#include <limits>
#include <utility>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cstdlib>
#include <ctime>
#include <random>
#include <chrono>

double w_ss = 1.08;
double w_ds = 0.34;
double w_c = std::pow(1.8, 2);
double w_b = 0.0001;
double w_ds_min = 10.09;
double initiation_dH = 0.2;
double initiation_dS = -5.7;
double terminal_AT_penalty_dH = 2.2;
double terminal_AT_penalty_dS = 6.9;
double R = 8.314;
double M = 1.0;
double k_plus = std::pow(10, 6);      
double tris = 40e-3;                                                  // [Tris]
double mg = 12.5e-3;                                                  // [Mg2+]
double staple_concentration = 100 * std::pow(10, -9);                 // [i]
double n = 0;                                                         // Coaxial stacking strength
double C = 2.8;                                                       // k+
double gamma_val = 2.5;                                               // Loop exponent
bool seq = true;                                                     // If false take average sequence free energy otherwise real sequence free energy
double T_min = 293.15;
double T_max = 373.15;
double constantT = 323.15;
double dTdt = 0.0166;
double dt = 1;
int n_traj = 1;
int threshold = 9;
double calculate_weight(const std::vector<int>& domain_info, int state) {
	int start_position = domain_info[1];
	int end_position = domain_info[2];

	// Calculate the weight based on the state
	if (state == 0) {
		return w_ss * (end_position - start_position + 1);
	}
	else {
		if (end_position - start_position + 1 > threshold) {
			double weight = w_ds * (end_position - start_position + 1);
			return std::pow(weight, 2);
		}
		else {
			return (w_ds_min);
		}
	}
}
std::vector<std::vector<double>> create_adjacency_matrix(const std::map<int, std::vector<std::vector<int>>>& staples, const std::map<int, std::vector<std::vector<int>>>& vertex, const std::map<int, std::vector<int>>& scafcross, const std::map<std::vector<int>, double>& dangling_ends) {
	int max_vertex = 0;

	// Find the maximum vertex
	for (const auto& vertex_entry : vertex) {
		for (const auto& domain_vertices : vertex_entry.second) {
			for (int vertex_num : domain_vertices) {
				if (vertex_num > max_vertex) {
					max_vertex = vertex_num;
				}
			}
		}
	}

	// Update maximum vertex considering dangling_ends map keys
	for (const auto& dangling_end_entry : dangling_ends) {
		for (int vertex_num : dangling_end_entry.first) {
			if (vertex_num > max_vertex) {
				max_vertex = vertex_num;
			}
		}
	}

	// Initialize the adjacency matrix
	std::vector<std::vector<double>> adjacency_matrix(max_vertex + 1, std::vector<double>(max_vertex + 1, 0.0));
	// Loop through the domains and calculate the weights
	for (const auto& staple : staples) {
		int staple_key = staple.first;
		const auto& domain_list = staple.second;

		for (size_t i = 0; i < domain_list.size(); i++) {
			const auto& domain = domain_list[i];
			int start_vertex = vertex.at(staple_key)[i][0];
			int end_vertex = vertex.at(staple_key)[i][1];
			double weight = calculate_weight(domain, 0);

			adjacency_matrix[start_vertex][end_vertex] = weight;
			adjacency_matrix[end_vertex][start_vertex] = weight;
		}
	}

	// Handle scaffold crossovers
	for (const auto& scafcross_entry : scafcross) {
		int helix_number = scafcross_entry.first;
		const auto& crossover_positions = scafcross_entry.second;

		for (int crossover_position : crossover_positions) {
			for (const auto& staple : staples) {
				const auto& domain_list = staple.second;

				for (size_t i = 0; i < domain_list.size(); i++) {
					const auto& domain = domain_list[i];

					// Find the two domains involved in the scaffold crossover
					if (domain[0] == helix_number && (domain[1] == crossover_position || domain[2] == crossover_position)) {
						int first_vertex;
						if (domain[1] == crossover_position) {
							first_vertex = vertex.at(staple.first)[i][0];
						}
						else {
							first_vertex = vertex.at(staple.first)[i][1];
						}

						int second_helix_number = helix_number + 1;
						for (const auto& second_staple : staples) {
							const auto& second_domain_list = second_staple.second;

							for (size_t j = 0; j < second_domain_list.size(); j++) {
								const auto& second_domain = second_domain_list[j];

								if (second_domain[0] == second_helix_number && (second_domain[1] == crossover_position || second_domain[2] == crossover_position)) {
									int second_vertex;
									if (second_domain[1] == crossover_position) {
										second_vertex = vertex.at(second_staple.first)[j][0];
									}
									else {
										second_vertex = vertex.at(second_staple.first)[j][1];
									}

									adjacency_matrix[first_vertex][second_vertex] += w_c;
									adjacency_matrix[second_vertex][first_vertex] += w_c;
								}
							}
						}
					}
				}
			}
		}
	}

	// Loop through the domains and calculate the weights
	// (Unchanged code)

	// Handle scaffold crossovers
	// (Unchanged code)

	// Loop through dangling_ends map and update the adjacency matrix
	for (const auto& dangling_end_entry : dangling_ends) {
		int start_vertex = dangling_end_entry.first[0];
		int end_vertex = dangling_end_entry.first[1];
		double weight = dangling_end_entry.second;

		adjacency_matrix[start_vertex][end_vertex] = weight;
		adjacency_matrix[end_vertex][start_vertex] = weight;
	}

	return adjacency_matrix;
}

struct CompareDist {
	bool operator()(const std::pair<int, double>& p1, const std::pair<int, double>& p2) {
		return p1.second > p2.second;
	}
};

double dijkstra_shortest_path(const std::vector<std::vector<double>>& adjacency_matrix, int start_vertex, int end_vertex) {
	int num_vertices = adjacency_matrix.size();
	std::vector<double> dist(num_vertices, std::numeric_limits<double>::max());
	std::vector<bool> visited(num_vertices, false);

	std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double>>, CompareDist> min_heap;
	dist[start_vertex] = 0;
	min_heap.push({ start_vertex, 0 });

	while (!min_heap.empty()) {
		int current_vertex = min_heap.top().first;
		min_heap.pop();

		if (visited[current_vertex]) {
			continue;
		}
		visited[current_vertex] = true;

		if (current_vertex == end_vertex) {
			break;
		}

		for (int neighbor = 0; neighbor < num_vertices; ++neighbor) {
			double edge_weight = adjacency_matrix[current_vertex][neighbor];
			if (edge_weight != 0.0 && !visited[neighbor]) {
				double new_dist = dist[current_vertex] + edge_weight;
				if (new_dist < dist[neighbor]) {
					dist[neighbor] = new_dist;
					min_heap.push({ neighbor, new_dist });
				}
			}
		}
	}

	if (visited[end_vertex]) {
		return dist[end_vertex];
	}
	else {
		return std::numeric_limits<double>::max();
	}
}

struct PairParams {
	double dH;
	double dS;
};

std::map<std::string, PairParams> nn_params = {
	{"AA", {-7.6, -21.3}},
	{"TT", {-7.6, -21.3}},
	{"AT", {-7.2, -20.4}},
	{"TA", {-7.2, -21.3}},
	{"CA", {-8.5, -22.7}},
	{"TG", {-8.5, -22.7}},
	{"GT", {-8.4, -22.4}},
	{"AC", {-8.4, -22.4}},
	{"CT", {-7.8, -21.0}},
	{"AG", {-7.8, -21.0}},
	{"GA", {-8.2, -22.2}},
	{"TC", {-8.2, -22.2}},
	{"CG", {-10.6, -27.2}},
	{"GC", {-9.8, -24.4}},
	{"GG", {-8.0, -19.9}},
	{"CC", {-8.0, -19.9}}

};

double calculate_santa_lucia_nn_free_energy(const std::vector<int>& dna_sequence, double temperature) {
	if (dna_sequence.size() < 2) {
		return 0.0;
	}

	std::vector<int> sequence = dna_sequence;

	std::vector<char> nucleotides_map = { 'A', 'T', 'C', 'G' };
	double total_dH = initiation_dH;
	double total_dS = initiation_dS;

	// Loop through the given DNA sequence and calculate the sum of dH and dS for each pair
	for (size_t i = 0; i < sequence.size() - 1; ++i) {
		char n1 = nucleotides_map[sequence[i] - 1];
		char n2 = nucleotides_map[sequence[i + 1] - 1];
		std::string nn_pair = { n1, n2 };
		total_dH += nn_params[nn_pair].dH;
		total_dS += nn_params[nn_pair].dS;
	}

	// Add terminal AT penalty if needed
	if (sequence.front() == 1 || sequence.front() == 2) {
		total_dH += terminal_AT_penalty_dH;
		total_dS += terminal_AT_penalty_dS;
	}
	if (sequence.back() == 1 || sequence.back() == 2) {
		total_dH += terminal_AT_penalty_dH;
		total_dS += terminal_AT_penalty_dS;
	}

	// Calculate free energy
	double dG = total_dH - temperature * total_dS / 1000; // in kcal/mol
	double energy_joules = dG * 4184; // Convert kcal/mol to J/mol

	// Calculate entropic penalty
	int num_phosphates = 2 * (sequence.size() - 1);
	double entropic_penalty = (0.368 * num_phosphates * std::log(0.5 * tris + 3.3 * std::sqrt(mg))) / 2;
	energy_joules -= temperature * entropic_penalty;
	return energy_joules;
}

int count_bound_neighbors(const std::map<int, std::vector<std::vector<std::vector<int>>>>& neighbors, const std::map<int, std::vector<int>>& state, int staple_id, int domain_idx) {
	int bound_neighbors = 0;

	// Get the neighbors of the given staple domain
	const std::vector<std::vector<std::vector<int>>>& domain_neighbors = neighbors.at(staple_id);

	// Loop through the neighbors of the given domain
	for (const std::vector<int>& neighbor_group : domain_neighbors[domain_idx]) {
		int neighbor_staple = neighbor_group[0];
		int neighbor_domain = neighbor_group[1];

		// Check if the neighboring domain is bound (state != 0)
		if (state.at(neighbor_staple)[neighbor_domain] != 0) {
			bound_neighbors++;
		}
	}
	return bound_neighbors;
}

double calculate_average_binding_energy(double temperature) {
	double sum_energy = 0;
	int count = 0;

	for (const auto& pair : nn_params) {
		double dH = pair.second.dH;
		double dS = pair.second.dS;
		double dG = dH - temperature * dS / 1000; // in kcal/mol
		double energy_joules = dG * 4184; // Convert kcal/mol to J/mol
		sum_energy += energy_joules;
		count++;
	}
	return sum_energy / count;
}

double calculate_unbinding_rate(double duplex_free_energy, int bound_neighbors, double temperature, double average_binding_energy) {


	// Calculate dGStack using the coaxial stacking strength n from the code
	double dGStack = n * bound_neighbors * average_binding_energy;
	// Calculate the rate of unbinding
	double rate_unbinding = k_plus * exp((duplex_free_energy + dGStack) / (R * temperature)) * M;

	return rate_unbinding;
}

double calculate_free_staple_binding_rate() {
	double rate = k_plus * staple_concentration;
	return rate;
}

double calculate_bound_staple_binding_rate(double shortest_path) {
	double E_r = shortest_path + w_c;
	double r_binding = k_plus * std::pow(C / E_r, static_cast<double>(gamma_val)) * M;
	return r_binding;
}

std::map<int, std::vector<double>> calculate_binding_energies(std::map<int, std::vector<std::vector<int>>>& sequence, std::map<int, std::vector<std::vector<int>>> stapPaths, double temperature, double average_binding_energy, bool seq) {
	std::map<int, std::vector<double>> binding_energy;

	// Loop through each staple using the map
	for (const auto& staple_pair : sequence) {
		int staple_id = staple_pair.first; // Get the staple_id from the map
		const auto& staple_domains = staple_pair.second; // Get the domains for the current staple

		std::vector<double> domain_binding_energies;

		// Loop through each domain of the current staple
		for (size_t i = 0; i < staple_domains.size(); i++) {
			double energy;

			if (seq) {
				energy = calculate_santa_lucia_nn_free_energy(staple_domains[i], temperature);
			}
			else {
				int domain_length = stapPaths[staple_id][i][2] - stapPaths[staple_id][i][1] + 1;
				energy = domain_length * average_binding_energy;
			}

			domain_binding_energies.push_back(energy);
		}

		// Add the binding energies for the domains of the current staple to the binding_energy map
		binding_energy[staple_id] = domain_binding_energies;
	}

	return binding_energy;
}

struct Event {
	double rate;
	int staple_id;
	std::vector<int> new_state;
	std::vector<double> modified_weights;
	std::vector<std::pair<int, int>> modified_edges;
};
/*
std::pair< std::vector<Event>, double> handle_events(std::map<int, std::vector<std::vector<int>>> stapPaths, std::map<int, std::vector<int>> con, std::map<int, std::vector<std::vector<int>>> vertex,
	std::map<int, std::vector<int>> state, std::map<int, std::vector<std::vector<std::vector<int>>>> neighbor, std::map<int, std::vector<double>> binding_energy,
	std::vector<std::vector<double>> adjacency_matrix, double T, double average_binding_energy) {
	std::vector<Event> events;
	double propensity = 0.0;
	for (const auto& staple_pair : stapPaths) {
		int staple_id = staple_pair.first; // Get the staple_id from the map

		std::vector<int>  staple_state = state[staple_id];
		std::vector<std::vector<int>> staple_path = stapPaths[staple_id];
		std::vector<std::vector<std::vector<int>>> staple_neighbors = neighbor[staple_id];
		std::vector<std::vector<int>> staple_vertex = vertex[staple_id];
		std::vector<double> staple_binding_energy = binding_energy[staple_id];
		std::vector<int> staple_connectivity = con[staple_id];
		// Loop through each domain of the current staple
		for (int domain_idx = 0; domain_idx < staple_state.size(); ++domain_idx) {

			int domain_state = staple_state[domain_idx];

			// If the domain_state is 0, it means it's a binding event
			if (domain_state == 0) {
				// Calculate the rate, new_state, modified_weights, and modified_edges
				// for the binding event

				// 1. Calculate the rate
				double rate = calculate_free_staple_binding_rate(); // Calculate the binding rate

				// 2. Update the state
				std::vector<int> new_state = staple_state;
				int max_state = *std::max_element(staple_state.begin(), staple_state.end());
				bool found_max = false;

				for (std::size_t i = 0; i < new_state.size(); ++i) {
					if (i == domain_idx) {
						// Check if the next domain or any domain to the right is equal to the max domain
						for (std::size_t j = i + 1; j < new_state.size(); ++j) {
							if (new_state[j] == max_state && new_state[j] != 0) {
								found_max = true;
								break;
							}
						}

						if (found_max) {
							new_state[i] = max_state;
						}
						else {
							new_state[i] = max_state + 1;
						}
					}
					else if (i > domain_idx && new_state[i] >= max_state && new_state[i] != 0) {
						// Increment the elements to the right of the current domain if they are greater than or equal to max_state and not equal to 0
						new_state[i]++;
					}
				}

				// 3. Calculate the modified weights
				std::vector<int> domain = staple_path[domain_idx];
				std::vector<double> modified_weights = { calculate_weight(domain, 1) };

				// 4. Determine the modified edges
				std::pair<int, int> edge_pair(staple_vertex.at(domain_idx)[0], staple_vertex.at(domain_idx)[1]);
				std::vector<std::pair<int, int>> modified_edges = { edge_pair };

				propensity += rate;

				// Create an Event with the calculated values
				Event event = { rate, staple_id, new_state, modified_weights, modified_edges };

				// Add the event to the events vector
				events.push_back(event);

				// Checking if the previous or next domain is already bound
				int prev_domain_idx = domain_idx - 1;
				int next_domain_idx = domain_idx + 1;

				bool prev_domain_bound = (prev_domain_idx >= 0) && (staple_state[prev_domain_idx] != 0);
				bool next_domain_bound = (next_domain_idx < staple_state.size()) && (staple_state[next_domain_idx] != 0);

				// Checking if the previous domain is already bound
				if (prev_domain_bound) {
					int case_connectivity = staple_connectivity[prev_domain_idx]; // Get the case from the connectivity vector using prev_domain_idx

					// Initialize new_state vector for the previous domain
					std::vector<int> new_state_previous = staple_state;
					new_state_previous[domain_idx] = staple_state[prev_domain_idx]; // Update the new_state_previous vector
					int start_vertex, end_vertex;

					// Initialize modified weights vector and modified edges vector for the previous domain
					std::vector<double> modified_weight_previous;
					std::vector<std::pair<int, int>> modified_edge_previous;

					// Determine the start and end vertices and weights based on the connectivity case

					if (case_connectivity == 0) {
						start_vertex = staple_vertex[prev_domain_idx][0];
						end_vertex = staple_vertex[domain_idx][0];

						// Modify the edges
						// 1. Include the edge for the crossover
						modified_edge_previous.push_back({ start_vertex, end_vertex });

						// 2. Include the edge for the current staple domain becoming bound
						std::pair<int, int> edge_pair(staple_vertex.at(domain_idx)[0], staple_vertex.at(domain_idx)[1]);
						modified_edge_previous.push_back(edge_pair);

						// Modify the weights
						// 1. Include the weight for the crossover
						modified_weight_previous.push_back(w_c);

						// 2. Include the weight for the current staple domain becoming bound
						double w_bound = calculate_weight(staple_path[domain_idx], 1);
						modified_weight_previous.push_back(w_bound);
					}
					else if (case_connectivity == 1) {
						start_vertex = staple_vertex[prev_domain_idx][1];
						end_vertex = staple_vertex[domain_idx][1];

						// Modify the edges
						// 1. Include the edge for the crossover
						modified_edge_previous.push_back({ start_vertex, end_vertex });

						// 2. Include the edge for the current staple domain becoming bound
						std::pair<int, int> edge_pair(staple_vertex.at(domain_idx)[0], staple_vertex.at(domain_idx)[1]);
						modified_edge_previous.push_back(edge_pair);

						// Modify the weights
						// 1. Include the weight for the crossover
						modified_weight_previous.push_back(w_c);

						// 2. Include the weight for the current staple domain becoming bound
						double w_bound = calculate_weight(staple_path[domain_idx], 1);
						modified_weight_previous.push_back(w_bound);

					}
					else if (case_connectivity == 2) {
						start_vertex = staple_vertex[prev_domain_idx][1];
						end_vertex = staple_vertex[domain_idx][0];

						// Modify the edges
						// 1. Include the edge for the bridge
						modified_edge_previous.push_back({ start_vertex, end_vertex });

						// 2. Include the edge for the current staple domain becoming bound
						std::pair<int, int> edge_pair(staple_vertex.at(domain_idx)[0], staple_vertex.at(domain_idx)[1]);
						modified_edge_previous.push_back(edge_pair);

						// Modify the weights
						// 1. Include the weight for the bridge
						modified_weight_previous.push_back(w_b);

						// 2. Include the weight for the current staple domain becoming bound
						double w_bound = calculate_weight(staple_path[domain_idx], 1);
						modified_weight_previous.push_back(w_bound);

					}
					else {
						start_vertex = staple_vertex[prev_domain_idx][0];
						end_vertex = staple_vertex[domain_idx][1];

						// Modify the edges
						// 1. Include the edge for the bridge
						modified_edge_previous.push_back({ start_vertex, end_vertex });

						// 2. Include the edge for the current staple domain becoming bound
						std::pair<int, int> edge_pair(staple_vertex.at(domain_idx)[0], staple_vertex.at(domain_idx)[1]);
						modified_edge_previous.push_back(edge_pair);

						// Modify the weights
						// 1. Include the weight for the crossover
						modified_weight_previous.push_back(w_b);

						// 2. Include the weight for the current staple domain becoming bound
						double w_bound = calculate_weight(staple_path[domain_idx], 1);
						modified_weight_previous.push_back(w_bound);
					}

					double distance = dijkstra_shortest_path(adjacency_matrix, start_vertex, end_vertex);
					double rate_previous = calculate_bound_staple_binding_rate(distance);

					propensity += rate_previous;

					Event event_previous = { rate_previous, staple_id, new_state_previous, modified_weight_previous, modified_edge_previous };

					// Add the event_previous to the events vector
					events.push_back(event_previous);
				}
				if (next_domain_bound) {
					int case_connectivity = staple_connectivity[domain_idx]; // Get the case from the connectivity vector using domain_idx
					std::vector<int> new_state_next = staple_state;
					new_state_next[domain_idx] = staple_state[next_domain_idx];  // Update the new_state_next vector

					int start_vertex, end_vertex;
					std::vector<double> modified_weight_next;
					std::vector<std::pair<int, int>> modified_edge_next;

					// Determine the start and end vertices and weights based on the connectivity case
					if (case_connectivity == 0) {
						start_vertex = staple_vertex[next_domain_idx][0];
						end_vertex = staple_vertex[domain_idx][0];

						// Modify the edges
						// 1. Include the edge for the crossover
						modified_edge_next.push_back({ start_vertex, end_vertex });

						// 2. Include the edge for the current staple domain becoming bound
						std::pair<int, int> edge_pair(staple_vertex.at(domain_idx)[0], staple_vertex.at(domain_idx)[1]);
						modified_edge_next.push_back(edge_pair);

						// Modify the weights
						// 1. Include the weight for the crossover
						modified_weight_next.push_back(w_c);

						// 2. Include the weight for the current staple domain becoming bound
						double w_bound = calculate_weight(staple_path[domain_idx], 1);
						modified_weight_next.push_back(w_bound);
					}
					else if (case_connectivity == 1) {
						start_vertex = staple_vertex[next_domain_idx][1];
						end_vertex = staple_vertex[domain_idx][1];

						// Modify the edges
						// 1. Include the edge for the crossover
						modified_edge_next.push_back({ start_vertex, end_vertex });

						// 2. Include the edge for the current staple domain becoming bound
						std::pair<int, int> edge_pair(staple_vertex.at(domain_idx)[0], staple_vertex.at(domain_idx)[1]);
						modified_edge_next.push_back(edge_pair);

						// Modify the weights
						// 1. Include the weight for the crossover
						modified_weight_next.push_back(w_c);

						// 2. Include the weight for the current staple domain becoming bound
						double w_bound = calculate_weight(staple_path[domain_idx], 1);
						modified_weight_next.push_back(w_bound);
					}
					else if (case_connectivity == 2) {
						start_vertex = staple_vertex[domain_idx][1];
						end_vertex = staple_vertex[next_domain_idx][0];

						// Modify the edges
						// 1. Include the edge for the bridge
						modified_edge_next.push_back({ start_vertex, end_vertex });

						// 2. Include the edge for the current staple domain becoming bound
						std::pair<int, int> edge_pair(staple_vertex.at(domain_idx)[0], staple_vertex.at(domain_idx)[1]);
						modified_edge_next.push_back(edge_pair);

						// Modify the weights
						// 1. Include the weight for the bridge
						modified_weight_next.push_back(w_b);

						// 2. Include the weight for the current staple domain becoming bound
						double w_bound = calculate_weight(staple_path[domain_idx], 1);
						modified_weight_next.push_back(w_bound);
					}
					else {
						start_vertex = staple_vertex[next_domain_idx][1];
						end_vertex = staple_vertex[domain_idx][0];

						// Modify the edges
						// 1. Include the edge for the bridge
						modified_edge_next.push_back({ start_vertex, end_vertex });

						// 2. Include the edge for the current staple domain becoming bound
						std::pair<int, int> edge_pair(staple_vertex.at(domain_idx)[0], staple_vertex.at(domain_idx)[1]);
						modified_edge_next.push_back(edge_pair);

						// Modify the weights
						// 1. Include the weight for the bridge
						modified_weight_next.push_back(w_b);

						// 2. Include the weight for the current staple domain becoming bound
						double w_bound = calculate_weight(staple_path[domain_idx], 1);
						modified_weight_next.push_back(w_bound);
					}

					double distance = dijkstra_shortest_path(adjacency_matrix, start_vertex, end_vertex);
					double rate_next = calculate_bound_staple_binding_rate(distance);

					propensity += rate_next;

					Event event_next = { rate_next, staple_id, new_state_next, modified_weight_next, modified_edge_next };

					// Add the event_next to the events vector
					events.push_back(event_next);
				}
			}
			else {

				int prev_domain_idx = domain_idx - 1;
				int next_domain_idx = domain_idx + 1;

				bool prev_domain_different = (prev_domain_idx >= 0) && (staple_state[prev_domain_idx] != staple_state[domain_idx]);
				bool next_domain_different = (next_domain_idx < staple_state.size()) && (staple_state[next_domain_idx] != staple_state[domain_idx]);
				if (domain_idx == 0 || prev_domain_different) {

					std::vector<int> new_state_previous = staple_state;
					new_state_previous[domain_idx] = 0;

					auto it = std::find(new_state_previous.begin(), new_state_previous.end(), staple_state[domain_idx]);

					if (it == new_state_previous.end()) {
						for (int i = 0; i < new_state_previous.size(); i++) {
							if (new_state_previous[i] > staple_state[domain_idx]) {
								new_state_previous[i]--;
							}
						}
					}
					std::vector<double> modified_weight_previous;
					std::vector<std::pair<int, int>> modified_edge_previous;

					if (domain_idx < staple_state.size() - 1 && staple_state[domain_idx] == staple_state[domain_idx + 1]) {
						// Remove the crossover
						int case_connectivity = staple_connectivity[domain_idx];
						// Handle the four cases
						if (case_connectivity == 0) {
							// Handle case 0
							int staple_id_neighbor_current = staple_neighbors[domain_idx][0][0];
							int domain_id_neighbor_current = staple_neighbors[domain_idx][0][1];
							int staple_id_neighbor_next = staple_neighbors[domain_idx + 1][0][0];
							int domain_id_neighbor_next = staple_neighbors[domain_idx + 1][0][1];

							if (staple_id_neighbor_current == staple_id_neighbor_next) {
								int current_domain_start = staple_path[domain_idx][1];
								int neighbor_domain_start = stapPaths[staple_id_neighbor_current][domain_id_neighbor_current][1];
								int connect_neighbor = std::min(domain_id_neighbor_current, domain_id_neighbor_next);
								int state_neighbor_current = state[staple_id_neighbor_current][domain_id_neighbor_current];
								int state_neighbor_next = state[staple_id_neighbor_next][domain_id_neighbor_next];

								if (!(con[staple_id_neighbor_current][connect_neighbor] == 1 && state_neighbor_current != 0 && state_neighbor_current == state_neighbor_next && current_domain_start > neighbor_domain_start)) {
									// Remove crossover
									int start_vertex = staple_vertex[domain_idx][0];
									int end_vertex = staple_vertex[domain_idx + 1][0];
									modified_edge_previous.push_back({ start_vertex, end_vertex });
									modified_weight_previous.push_back(0);
								}
							}
						}
						else if (case_connectivity == 1) {
							// Handle case 1
							int staple_id_neighbor_current = staple_neighbors[domain_idx].back()[0];
							int domain_id_neighbor_current = staple_neighbors[domain_idx].back()[1];
							int staple_id_neighbor_next = staple_neighbors[domain_idx + 1].back()[0];
							int domain_id_neighbor_next = staple_neighbors[domain_idx + 1].back()[1];

							if (staple_id_neighbor_current == staple_id_neighbor_next) {
								int current_domain_start = staple_path[domain_idx][1];
								int neighbor_domain_start = stapPaths[staple_id_neighbor_current][domain_id_neighbor_current][1];
								int connect_neighbor = std::min(domain_id_neighbor_current, domain_id_neighbor_next);
								int state_neighbor_current = state[staple_id_neighbor_current][domain_id_neighbor_current];
								int state_neighbor_next = state[staple_id_neighbor_next][domain_id_neighbor_next];

								if (!(con[staple_id_neighbor_current][connect_neighbor] == 0 && state_neighbor_current != 0 && state_neighbor_current == state_neighbor_next && current_domain_start < neighbor_domain_start)) {
									// Remove crossover
									int start_vertex = staple_vertex[domain_idx][1];
									int end_vertex = staple_vertex[domain_idx + 1][1];
									modified_edge_previous.push_back({ start_vertex, end_vertex });
									modified_weight_previous.push_back(0);
								}
							}
						}
						else if (case_connectivity == 2) {
							int start_vertex = staple_vertex[domain_idx][1];
							int end_vertex = staple_vertex[domain_idx + 1][0];
							modified_edge_previous.push_back({ start_vertex, end_vertex });
							modified_weight_previous.push_back(0);
						}
						else {
							int start_vertex = staple_vertex[domain_idx + 1][1];
							int end_vertex = staple_vertex[domain_idx][0];
							modified_edge_previous.push_back({ start_vertex, end_vertex });
							modified_weight_previous.push_back(0);
						}

					}

					std::pair<int, int> edge_pair(staple_vertex.at(domain_idx)[0], staple_vertex.at(domain_idx)[1]);
					modified_edge_previous.push_back(edge_pair);

					double w_bound_previous = calculate_weight(staple_path[domain_idx], 0);
					modified_weight_previous.push_back(w_bound_previous);

					int bound_neighbors = count_bound_neighbors(neighbor, state, staple_id, domain_idx);
					double rate_previous = calculate_unbinding_rate(staple_binding_energy[domain_idx], bound_neighbors, T, average_binding_energy);

					propensity += rate_previous;

					Event event_previous = { rate_previous, staple_id, new_state_previous, modified_weight_previous, modified_edge_previous };

					// Add the event to the events vector
					events.push_back(event_previous);

				}
				if ((domain_idx == staple_state.size() - 1 || next_domain_different) && !(domain_idx == 0 || prev_domain_different)) {
					// Initialize new_state_next, modified_weight_next, and modified_edge_next
					std::vector<int> new_state_next = staple_state;
					new_state_next[domain_idx] = 0;

					auto it = std::find(new_state_next.begin(), new_state_next.end(), staple_state[domain_idx]);

					if (it == new_state_next.end()) {
						for (int i = 0; i < new_state_next.size(); i++) {
							if (new_state_next[i] > staple_state[domain_idx]) {
								new_state_next[i]--;
							}
						}
					}

					std::vector<double> modified_weight_next;
					std::vector<std::pair<int, int>> modified_edge_next;

					if (domain_idx > 0 && staple_state[domain_idx] == staple_state[domain_idx - 1]) {
						// The previous domain is bound to the current domain
						int case_connectivity = staple_connectivity[domain_idx - 1];
						if (case_connectivity == 0) {
							// Handle case 0
							int staple_id_neighbor_prev = staple_neighbors[domain_idx - 1][0][0];
							int domain_id_neighbor_prev = staple_neighbors[domain_idx - 1][0][1];
							int staple_id_neighbor_current = staple_neighbors[domain_idx][0][0];
							int domain_id_neighbor_current = staple_neighbors[domain_idx][0][1];

							if (staple_id_neighbor_prev == staple_id_neighbor_current) {
								int current_domain_start = staple_path[domain_idx][1];
								int neighbor_domain_start = stapPaths[staple_id_neighbor_current][domain_id_neighbor_current][1];
								int connect_neighbor = std::min(domain_id_neighbor_prev, domain_id_neighbor_current);
								int state_neighbor_prev = state[staple_id_neighbor_prev][domain_id_neighbor_prev];
								int state_neighbor_current = state[staple_id_neighbor_current][domain_id_neighbor_current];

								if (!(con[staple_id_neighbor_prev][connect_neighbor] == 1 && state_neighbor_prev != 0 && state_neighbor_prev == state_neighbor_current && current_domain_start > neighbor_domain_start)) {
									// Remove crossover
									int start_vertex = staple_vertex[domain_idx - 1][0];
									int end_vertex = staple_vertex[domain_idx][0];
									modified_edge_next.push_back({ start_vertex, end_vertex });
									modified_weight_next.push_back(0);
								}
							}
						}
						else if (case_connectivity == 1) {
							// Handle case 1
							int staple_id_neighbor_prev = staple_neighbors[domain_idx - 1].back()[0];
							int domain_id_neighbor_prev = staple_neighbors[domain_idx - 1].back()[1];
							int staple_id_neighbor_current = staple_neighbors[domain_idx].back()[0];
							int domain_id_neighbor_current = staple_neighbors[domain_idx].back()[1];

							if (staple_id_neighbor_prev == staple_id_neighbor_current) {
								int current_domain_start = staple_path[domain_idx][1];
								int neighbor_domain_start = stapPaths[staple_id_neighbor_current][domain_id_neighbor_current][1];
								int connect_neighbor = std::min(domain_id_neighbor_prev, domain_id_neighbor_current);
								int state_neighbor_prev = state[staple_id_neighbor_prev][domain_id_neighbor_prev];
								int state_neighbor_current = state[staple_id_neighbor_current][domain_id_neighbor_current];

								if (!(con[staple_id_neighbor_prev][connect_neighbor] == 0 && state_neighbor_prev != 0 && state_neighbor_prev == state_neighbor_current && current_domain_start < neighbor_domain_start)) {
									// Remove crossover
									int start_vertex = staple_vertex[domain_idx - 1][1];
									int end_vertex = staple_vertex[domain_idx][1];
									modified_edge_next.push_back({ start_vertex, end_vertex });
									modified_weight_next.push_back(0);
								}
							}
						}
						else if (case_connectivity == 2) {
							// Handle case 2
							int start_vertex = staple_vertex[domain_idx - 1][1];
							int end_vertex = staple_vertex[domain_idx][0];
							modified_edge_next.push_back({ start_vertex, end_vertex });
							modified_weight_next.push_back(0);
						}
						else {
							// Handle case 3
							int start_vertex = staple_vertex[domain_idx][1];
							int end_vertex = staple_vertex[domain_idx - 1][0];
							modified_edge_next.push_back({ start_vertex, end_vertex });
							modified_weight_next.push_back(0);
						}
					}
					std::pair<int, int> edge_pair_next(staple_vertex.at(domain_idx)[0], staple_vertex.at(domain_idx)[1]);
					modified_edge_next.push_back(edge_pair_next);

					double w_bound_next = calculate_weight(staple_path[domain_idx], 0);
					modified_weight_next.push_back(w_bound_next);

					int bound_neighbors_next = count_bound_neighbors(neighbor, state, staple_id, domain_idx);
					double rate_next = calculate_unbinding_rate(staple_binding_energy[domain_idx], bound_neighbors_next, T, average_binding_energy);

					propensity += rate_next;

					Event event_next = { rate_next, staple_id, new_state_next, modified_weight_next, modified_edge_next };

					events.push_back(event_next);
				}
			}
		}
	}
	// Iterate through the staple_state vector to identify binding events
	return std::pair<std::vector<Event>, double>(events, propensity);
}
*/
void print_events(const std::vector<Event>& events) {
	for (const auto& event : events) {
		std::cout << "Event:\n";
		std::cout << "\trate: " << event.rate << "\n";
		std::cout << "\tstaple_id: " << event.staple_id << "\n";

		std::cout << "\tnew_state: ";
		for (const auto& state : event.new_state) {
			std::cout << state << " ";
		}
		std::cout << "\n";

		std::cout << "\tmodified_weights: ";
		for (const auto& weight : event.modified_weights) {
			std::cout << weight << " ";
		}
		std::cout << "\n";

		std::cout << "\tmodified_edges: ";
		for (const auto& edge : event.modified_edges) {
			std::cout << "(" << edge.first << ", " << edge.second << ") ";
		}
		std::cout << "\n";
	}
}
std::pair< std::vector<Event>, double> handle_events(std::map<int, std::vector<std::vector<int>>> stapPaths, std::map<int, std::vector<int>> con, std::map<int, std::vector<std::vector<int>>> vertex,
	std::map<int, std::vector<int>> state, std::map<int, std::vector<std::vector<std::vector<int>>>> neighbor, std::map<int, std::vector<double>> binding_energy,
	std::vector<std::vector<double>> adjacency_matrix, double T, double average_binding_energy) {
	std::vector<Event> events;
	double propensity = 0.0;
	for (const auto& staple_pair : stapPaths) {
		int staple_id = staple_pair.first; // Get the staple_id from the map

		std::vector<int>  staple_state = state[staple_id];
		std::vector<std::vector<int>> staple_path = stapPaths[staple_id];
		std::vector<std::vector<std::vector<int>>> staple_neighbors = neighbor[staple_id];
		std::vector<std::vector<int>> staple_vertex = vertex[staple_id];
		std::vector<double> staple_binding_energy = binding_energy[staple_id];
		std::vector<int> staple_connectivity = con[staple_id];
		// Loop through each domain of the current staple
		for (int domain_idx = 0; domain_idx < staple_state.size(); ++domain_idx) {

			int domain_state = staple_state[domain_idx];

			// If the domain_state is 0, it means it's a binding event
			if (domain_state == 0) {
				// Calculate the rate, new_state, modified_weights, and modified_edges
				// for the binding event

				// 1. Calculate the rate
				double rate = calculate_free_staple_binding_rate(); // Calculate the binding rate

				// 2. Update the state
				std::vector<int> new_state = staple_state;
				int max_state = *std::max_element(staple_state.begin(), staple_state.end());
				bool found_max = false;

				for (std::size_t i = 0; i < new_state.size(); ++i) {
					if (i == domain_idx) {
						// Check if the next domain or any domain to the right is equal to the max domain
						for (std::size_t j = i + 1; j < new_state.size(); ++j) {
							if (new_state[j] == max_state && new_state[j] != 0) {
								found_max = true;
								break;
							}
						}

						if (found_max) {
							new_state[i] = max_state;
						}
						else {
							new_state[i] = max_state + 1;
						}
					}
					else if (i > domain_idx && new_state[i] >= max_state && new_state[i] != 0) {
						// Increment the elements to the right of the current domain if they are greater than or equal to max_state and not equal to 0
						new_state[i]++;
					}
				}

				// 3. Calculate the modified weights
				std::vector<int> domain = staple_path[domain_idx];
				std::vector<double> modified_weights = { calculate_weight(domain, 1) };

				// 4. Determine the modified edges
				std::pair<int, int> edge_pair(staple_vertex.at(domain_idx)[0], staple_vertex.at(domain_idx)[1]);
				std::vector<std::pair<int, int>> modified_edges = { edge_pair };

				propensity += rate;

				// Create an Event with the calculated values
				Event event = { rate, staple_id, new_state, modified_weights, modified_edges };

				// Add the event to the events vector
				events.push_back(event);

				// Checking if the previous or next domain is already bound
				int prev_domain_idx = domain_idx - 1;
				int next_domain_idx = domain_idx + 1;

				bool prev_domain_bound = (prev_domain_idx >= 0) && (staple_state[prev_domain_idx] != 0);
				bool next_domain_bound = (next_domain_idx < staple_state.size()) && (staple_state[next_domain_idx] != 0);

				// Checking if the previous domain is already bound
				if (prev_domain_bound) {
					int case_connectivity = staple_connectivity[prev_domain_idx]; // Get the case from the connectivity vector using prev_domain_idx

					// Initialize new_state vector for the previous domain
					std::vector<int> new_state_previous = staple_state;
					new_state_previous[domain_idx] = staple_state[prev_domain_idx]; // Update the new_state_previous vector
					int start_vertex, end_vertex;

					// Initialize modified weights vector and modified edges vector for the previous domain
					std::vector<double> modified_weight_previous;
					std::vector<std::pair<int, int>> modified_edge_previous;

					// Determine the start and end vertices and weights based on the connectivity case 

					if (case_connectivity == 0) {
						start_vertex = staple_vertex[prev_domain_idx][0];
						end_vertex = staple_vertex[domain_idx][0];

						// Modify the edges
						// 1. Include the edge for the crossover
						modified_edge_previous.push_back({ start_vertex, end_vertex });

						// 2. Include the edge for the current staple domain becoming bound
						std::pair<int, int> edge_pair(staple_vertex.at(domain_idx)[0], staple_vertex.at(domain_idx)[1]);
						modified_edge_previous.push_back(edge_pair);

						// Modify the weights
						// 1. Include the weight for the crossover
						modified_weight_previous.push_back(w_c);

						// 2. Include the weight for the current staple domain becoming bound
						double w_bound = calculate_weight(staple_path[domain_idx], 1);
						modified_weight_previous.push_back(w_bound);
					}
					else if (case_connectivity == 1) {
						start_vertex = staple_vertex[prev_domain_idx][1];
						end_vertex = staple_vertex[domain_idx][1];

						// Modify the edges
						// 1. Include the edge for the crossover
						modified_edge_previous.push_back({ start_vertex, end_vertex });

						// 2. Include the edge for the current staple domain becoming bound
						std::pair<int, int> edge_pair(staple_vertex.at(domain_idx)[0], staple_vertex.at(domain_idx)[1]);
						modified_edge_previous.push_back(edge_pair);

						// Modify the weights
						// 1. Include the weight for the crossover
						modified_weight_previous.push_back(w_c);

						// 2. Include the weight for the current staple domain becoming bound
						double w_bound = calculate_weight(staple_path[domain_idx], 1);
						modified_weight_previous.push_back(w_bound);

					}
					else if (case_connectivity == 2) {
						start_vertex = staple_vertex[prev_domain_idx][1];
						end_vertex = staple_vertex[domain_idx][0];

						// Modify the edges
						// 1. Include the edge for the bridge
						modified_edge_previous.push_back({ start_vertex, end_vertex });

						// 2. Include the edge for the current staple domain becoming bound
						std::pair<int, int> edge_pair(staple_vertex.at(domain_idx)[0], staple_vertex.at(domain_idx)[1]);
						modified_edge_previous.push_back(edge_pair);

						// Modify the weights
						// 1. Include the weight for the bridge
						modified_weight_previous.push_back(w_b);

						// 2. Include the weight for the current staple domain becoming bound
						double w_bound = calculate_weight(staple_path[domain_idx], 1);
						modified_weight_previous.push_back(w_bound);

					}
					else {
						start_vertex = staple_vertex[prev_domain_idx][0];
						end_vertex = staple_vertex[domain_idx][1];

						// Modify the edges
						// 1. Include the edge for the bridge
						modified_edge_previous.push_back({ start_vertex, end_vertex });

						// 2. Include the edge for the current staple domain becoming bound
						std::pair<int, int> edge_pair(staple_vertex.at(domain_idx)[0], staple_vertex.at(domain_idx)[1]);
						modified_edge_previous.push_back(edge_pair);

						// Modify the weights
						// 1. Include the weight for the crossover
						modified_weight_previous.push_back(w_b);

						// 2. Include the weight for the current staple domain becoming bound
						double w_bound = calculate_weight(staple_path[domain_idx], 1);
						modified_weight_previous.push_back(w_bound);
					}

					double distance = dijkstra_shortest_path(adjacency_matrix, start_vertex, end_vertex);
					double rate_previous = calculate_bound_staple_binding_rate(distance);

					propensity += rate_previous;

					Event event_previous = { rate_previous, staple_id, new_state_previous, modified_weight_previous, modified_edge_previous };

					// Add the event_previous to the events vector
					events.push_back(event_previous);
				}
				if (next_domain_bound) {
					int case_connectivity = staple_connectivity[domain_idx]; // Get the case from the connectivity vector using domain_idx
					std::vector<int> new_state_next = staple_state;
					new_state_next[domain_idx] = staple_state[next_domain_idx];  // Update the new_state_next vector

					int start_vertex, end_vertex;
					std::vector<double> modified_weight_next;
					std::vector<std::pair<int, int>> modified_edge_next;

					// Determine the start and end vertices and weights based on the connectivity case
					if (case_connectivity == 0) {
						start_vertex = staple_vertex[next_domain_idx][0];
						end_vertex = staple_vertex[domain_idx][0];

						// Modify the edges
						// 1. Include the edge for the crossover
						modified_edge_next.push_back({ start_vertex, end_vertex });

						// 2. Include the edge for the current staple domain becoming bound
						std::pair<int, int> edge_pair(staple_vertex.at(domain_idx)[0], staple_vertex.at(domain_idx)[1]);
						modified_edge_next.push_back(edge_pair);

						// Modify the weights
						// 1. Include the weight for the crossover
						modified_weight_next.push_back(w_c);

						// 2. Include the weight for the current staple domain becoming bound
						double w_bound = calculate_weight(staple_path[domain_idx], 1);
						modified_weight_next.push_back(w_bound);
					}
					else if (case_connectivity == 1) {
						start_vertex = staple_vertex[next_domain_idx][1];
						end_vertex = staple_vertex[domain_idx][1];

						// Modify the edges
						// 1. Include the edge for the crossover
						modified_edge_next.push_back({ start_vertex, end_vertex });

						// 2. Include the edge for the current staple domain becoming bound
						std::pair<int, int> edge_pair(staple_vertex.at(domain_idx)[0], staple_vertex.at(domain_idx)[1]);
						modified_edge_next.push_back(edge_pair);

						// Modify the weights
						// 1. Include the weight for the crossover
						modified_weight_next.push_back(w_c);

						// 2. Include the weight for the current staple domain becoming bound
						double w_bound = calculate_weight(staple_path[domain_idx], 1);
						modified_weight_next.push_back(w_bound);
					}
					else if (case_connectivity == 2) {
						start_vertex = staple_vertex[domain_idx][1];
						end_vertex = staple_vertex[next_domain_idx][0];

						// Modify the edges
						// 1. Include the edge for the bridge
						modified_edge_next.push_back({ start_vertex, end_vertex });

						// 2. Include the edge for the current staple domain becoming bound
						std::pair<int, int> edge_pair(staple_vertex.at(domain_idx)[0], staple_vertex.at(domain_idx)[1]);
						modified_edge_next.push_back(edge_pair);

						// Modify the weights
						// 1. Include the weight for the bridge
						modified_weight_next.push_back(w_b);

						// 2. Include the weight for the current staple domain becoming bound
						double w_bound = calculate_weight(staple_path[domain_idx], 1);
						modified_weight_next.push_back(w_bound);
					}
					else {
						start_vertex = staple_vertex[next_domain_idx][1];
						end_vertex = staple_vertex[domain_idx][0];

						// Modify the edges
						// 1. Include the edge for the bridge
						modified_edge_next.push_back({ start_vertex, end_vertex });

						// 2. Include the edge for the current staple domain becoming bound
						std::pair<int, int> edge_pair(staple_vertex.at(domain_idx)[0], staple_vertex.at(domain_idx)[1]);
						modified_edge_next.push_back(edge_pair);

						// Modify the weights
						// 1. Include the weight for the bridge
						modified_weight_next.push_back(w_b);

						// 2. Include the weight for the current staple domain becoming bound
						double w_bound = calculate_weight(staple_path[domain_idx], 1);
						modified_weight_next.push_back(w_bound);
					}

					double distance = dijkstra_shortest_path(adjacency_matrix, start_vertex, end_vertex);
					double rate_next = calculate_bound_staple_binding_rate(distance);

					propensity += rate_next;

					Event event_next = { rate_next, staple_id, new_state_next, modified_weight_next, modified_edge_next };

					// Add the event_next to the events vector
					events.push_back(event_next);
				}
			}
			else {

				int prev_domain_idx = domain_idx - 1;
				int next_domain_idx = domain_idx + 1;
				bool prev_domain_different = (prev_domain_idx >= 0) && (staple_state[prev_domain_idx] != staple_state[domain_idx]);
				bool next_domain_different = (next_domain_idx < staple_state.size()) && (staple_state[next_domain_idx] != staple_state[domain_idx]);
				if (domain_idx == 0 || prev_domain_different) {
					std::vector<int> new_state_previous = staple_state;
					new_state_previous[domain_idx] = 0;

					auto it = std::find(new_state_previous.begin(), new_state_previous.end(), staple_state[domain_idx]);

					if (it == new_state_previous.end()) {
						for (int i = 0; i < new_state_previous.size(); i++) {
							if (new_state_previous[i] > staple_state[domain_idx]) {
								new_state_previous[i]--;
							}
						}
					}
					std::vector<double> modified_weight_previous;
					std::vector<std::pair<int, int>> modified_edge_previous;

					if (domain_idx < staple_state.size() - 1 && staple_state[domain_idx] == staple_state[domain_idx + 1]) {
						// Remove the crossover
						int case_connectivity = staple_connectivity[domain_idx];
						// Handle the four cases
						if (case_connectivity == 0) {
							// Handle case 0
							if (staple_neighbors[domain_idx].size() > 0) {
								int staple_id_neighbor_current = staple_neighbors[domain_idx][0][0];
								int domain_id_neighbor_current = staple_neighbors[domain_idx][0][1];
								int staple_id_neighbor_next = staple_neighbors[domain_idx + 1][0][0];
								int domain_id_neighbor_next = staple_neighbors[domain_idx + 1][0][1];
								if (staple_id_neighbor_current == staple_id_neighbor_next) {
									int current_domain_start = staple_path[domain_idx][1];
									int neighbor_domain_start = stapPaths[staple_id_neighbor_current][domain_id_neighbor_current][1];
									int connect_neighbor = std::min(domain_id_neighbor_current, domain_id_neighbor_next);
									int state_neighbor_current = state[staple_id_neighbor_current][domain_id_neighbor_current];
									int state_neighbor_next = state[staple_id_neighbor_next][domain_id_neighbor_next];

									if (!(con[staple_id_neighbor_current][connect_neighbor] == 1 && state_neighbor_current != 0 && state_neighbor_current == state_neighbor_next && current_domain_start > neighbor_domain_start)) {
										// Remove crossover
										int start_vertex = staple_vertex[domain_idx][0];
										int end_vertex = staple_vertex[domain_idx + 1][0];
										modified_edge_previous.push_back({ start_vertex, end_vertex });
										modified_weight_previous.push_back(0);
									}
								}

							}
							else {
								int start_vertex = staple_vertex[domain_idx][0];
								int end_vertex = staple_vertex[domain_idx + 1][0];
								modified_edge_previous.push_back({ start_vertex, end_vertex });
								modified_weight_previous.push_back(0);
							}
						}
						else if (case_connectivity == 1) {
							// Handle case 1
							if (staple_neighbors[domain_idx].size() > 0) {
								int staple_id_neighbor_current = staple_neighbors[domain_idx].back()[0];
								int domain_id_neighbor_current = staple_neighbors[domain_idx].back()[1];
								int staple_id_neighbor_next = staple_neighbors[domain_idx + 1].back()[0];
								int domain_id_neighbor_next = staple_neighbors[domain_idx + 1].back()[1];

								if (staple_id_neighbor_current == staple_id_neighbor_next) {
									int current_domain_start = staple_path[domain_idx][1];
									int neighbor_domain_start = stapPaths[staple_id_neighbor_current][domain_id_neighbor_current][1];
									int connect_neighbor = std::min(domain_id_neighbor_current, domain_id_neighbor_next);
									int state_neighbor_current = state[staple_id_neighbor_current][domain_id_neighbor_current];
									int state_neighbor_next = state[staple_id_neighbor_next][domain_id_neighbor_next];

									if (!(con[staple_id_neighbor_current][connect_neighbor] == 0 && state_neighbor_current != 0 && state_neighbor_current == state_neighbor_next && current_domain_start < neighbor_domain_start)) {
										// Remove crossover
										int start_vertex = staple_vertex[domain_idx][1];
										int end_vertex = staple_vertex[domain_idx + 1][1];
										modified_edge_previous.push_back({ start_vertex, end_vertex });
										modified_weight_previous.push_back(0);
									}
								}
							}
							else {
								int start_vertex = staple_vertex[domain_idx][1];
								int end_vertex = staple_vertex[domain_idx + 1][1];
								modified_edge_previous.push_back({ start_vertex, end_vertex });
								modified_weight_previous.push_back(0);
							}
						}
						else if (case_connectivity == 2) {
							int start_vertex = staple_vertex[domain_idx][1];
							int end_vertex = staple_vertex[domain_idx + 1][0];
							modified_edge_previous.push_back({ start_vertex, end_vertex });
							modified_weight_previous.push_back(0);
						}
						else {
							int start_vertex = staple_vertex[domain_idx + 1][1];
							int end_vertex = staple_vertex[domain_idx][0];
							modified_edge_previous.push_back({ start_vertex, end_vertex });
							modified_weight_previous.push_back(0);
						}

					}

					std::pair<int, int> edge_pair(staple_vertex.at(domain_idx)[0], staple_vertex.at(domain_idx)[1]);
					modified_edge_previous.push_back(edge_pair);

					double w_bound_previous = calculate_weight(staple_path[domain_idx], 0);
					modified_weight_previous.push_back(w_bound_previous);

					int bound_neighbors = count_bound_neighbors(neighbor, state, staple_id, domain_idx);
					double rate_previous = calculate_unbinding_rate(staple_binding_energy[domain_idx], bound_neighbors, T, average_binding_energy);

					propensity += rate_previous;

					Event event_previous = { rate_previous, staple_id, new_state_previous, modified_weight_previous, modified_edge_previous };

					// Add the event to the events vector
					events.push_back(event_previous);

				}
				if ((domain_idx == staple_state.size() - 1 || next_domain_different) && !(domain_idx == 0 || prev_domain_different)) {
					// Initialize new_state_next, modified_weight_next, and modified_edge_next
					std::vector<int> new_state_next = staple_state;
					new_state_next[domain_idx] = 0;

					auto it = std::find(new_state_next.begin(), new_state_next.end(), staple_state[domain_idx]);

					if (it == new_state_next.end()) {
						for (int i = 0; i < new_state_next.size(); i++) {
							if (new_state_next[i] > staple_state[domain_idx]) {
								new_state_next[i]--;
							}
						}
					}

					std::vector<double> modified_weight_next;
					std::vector<std::pair<int, int>> modified_edge_next;

					if (domain_idx > 0 && staple_state[domain_idx] == staple_state[domain_idx - 1]) {
						// The previous domain is bound to the current domain
						int case_connectivity = staple_connectivity[domain_idx - 1];
						if (case_connectivity == 0) {
							// Handle case 0
							if (staple_neighbors[domain_idx].size() > 0) {
								int staple_id_neighbor_prev = staple_neighbors[domain_idx - 1][0][0];
								int domain_id_neighbor_prev = staple_neighbors[domain_idx - 1][0][1];
								int staple_id_neighbor_current = staple_neighbors[domain_idx][0][0];
								int domain_id_neighbor_current = staple_neighbors[domain_idx][0][1];

								if (staple_id_neighbor_prev == staple_id_neighbor_current) {
									int current_domain_start = staple_path[domain_idx][1];
									int neighbor_domain_start = stapPaths[staple_id_neighbor_current][domain_id_neighbor_current][1];
									int connect_neighbor = std::min(domain_id_neighbor_prev, domain_id_neighbor_current);
									int state_neighbor_prev = state[staple_id_neighbor_prev][domain_id_neighbor_prev];
									int state_neighbor_current = state[staple_id_neighbor_current][domain_id_neighbor_current];

									if (!(con[staple_id_neighbor_prev][connect_neighbor] == 1 && state_neighbor_prev != 0 && state_neighbor_prev == state_neighbor_current && current_domain_start > neighbor_domain_start)) {
										// Remove crossover
										int start_vertex = staple_vertex[domain_idx - 1][0];
										int end_vertex = staple_vertex[domain_idx][0];
										modified_edge_next.push_back({ start_vertex, end_vertex });
										modified_weight_next.push_back(0);
									}
								}
							}
							else {
								int start_vertex = staple_vertex[domain_idx - 1][0];
								int end_vertex = staple_vertex[domain_idx][0];
								modified_edge_next.push_back({ start_vertex, end_vertex });
								modified_weight_next.push_back(0);
							}
						}
						else if (case_connectivity == 1) {
							// Handle case 1
							if (staple_neighbors[domain_idx].size() > 0) {
								int staple_id_neighbor_prev = staple_neighbors[domain_idx - 1].back()[0];
								int domain_id_neighbor_prev = staple_neighbors[domain_idx - 1].back()[1];
								int staple_id_neighbor_current = staple_neighbors[domain_idx].back()[0];
								int domain_id_neighbor_current = staple_neighbors[domain_idx].back()[1];

								if (staple_id_neighbor_prev == staple_id_neighbor_current) {
									int current_domain_start = staple_path[domain_idx][1];
									int neighbor_domain_start = stapPaths[staple_id_neighbor_current][domain_id_neighbor_current][1];
									int connect_neighbor = std::min(domain_id_neighbor_prev, domain_id_neighbor_current);
									int state_neighbor_prev = state[staple_id_neighbor_prev][domain_id_neighbor_prev];
									int state_neighbor_current = state[staple_id_neighbor_current][domain_id_neighbor_current];

									if (!(con[staple_id_neighbor_prev][connect_neighbor] == 0 && state_neighbor_prev != 0 && state_neighbor_prev == state_neighbor_current && current_domain_start < neighbor_domain_start)) {
										// Remove crossover
										int start_vertex = staple_vertex[domain_idx - 1][1];
										int end_vertex = staple_vertex[domain_idx][1];
										modified_edge_next.push_back({ start_vertex, end_vertex });
										modified_weight_next.push_back(0);
									}
								}
							}
							else {
								int start_vertex = staple_vertex[domain_idx - 1][1];
								int end_vertex = staple_vertex[domain_idx][1];
								modified_edge_next.push_back({ start_vertex, end_vertex });
								modified_weight_next.push_back(0);
							}
						}
						else if (case_connectivity == 2) {
							// Handle case 2
							int start_vertex = staple_vertex[domain_idx - 1][1];
							int end_vertex = staple_vertex[domain_idx][0];
							modified_edge_next.push_back({ start_vertex, end_vertex });
							modified_weight_next.push_back(0);
						}
						else {
							// Handle case 3
							int start_vertex = staple_vertex[domain_idx][1];
							int end_vertex = staple_vertex[domain_idx - 1][0];
							modified_edge_next.push_back({ start_vertex, end_vertex });
							modified_weight_next.push_back(0);
						}
					}
					std::pair<int, int> edge_pair_next(staple_vertex.at(domain_idx)[0], staple_vertex.at(domain_idx)[1]);
					modified_edge_next.push_back(edge_pair_next);

					double w_bound_next = calculate_weight(staple_path[domain_idx], 0);
					modified_weight_next.push_back(w_bound_next);

					int bound_neighbors_next = count_bound_neighbors(neighbor, state, staple_id, domain_idx);
					double rate_next = calculate_unbinding_rate(staple_binding_energy[domain_idx], bound_neighbors_next, T, average_binding_energy);

					propensity += rate_next;

					Event event_next = { rate_next, staple_id, new_state_next, modified_weight_next, modified_edge_next };

					events.push_back(event_next);
				}
			}
		}
	}
	// Iterate through the staple_state vector to identify binding events
	//print_events(events);
	return std::pair<std::vector<Event>, double>(events, propensity);
}

void apply_changes(const Event& event, std::map<int, std::vector<int>>& state, std::vector<std::vector<double>>& adjacency_matrix) {
	// Update the state of the staple
	state[event.staple_id] = event.new_state;
	//std::cout << event.new_state[0] << "--" << event.new_state[1] << "--" << event.new_state[2] << std::endl;
	// Apply the modified weights
	for (size_t i = 0; i < event.modified_edges.size(); i++) {
		const auto& edge = event.modified_edges[i];
		double weight = event.modified_weights[i];

		// Update the adjacency matrix for both directions (undirected weighted graph)
		adjacency_matrix[edge.first][edge.second] = weight;
		adjacency_matrix[edge.second][edge.first] = weight;
	}
}

std::vector<double> generate_temp_ramp(double high_temperature, double low_temperature, double cooling_rate) {
	std::vector<double> temperature_profile;

	// Calculate the time it takes to cool down
	int cooling_time = static_cast<int>((high_temperature - low_temperature) / cooling_rate);

	// Cool down
	for (int i = 0; i <= cooling_time; i++) {
		double current_temperature = high_temperature - i * cooling_rate;
		temperature_profile.push_back(current_temperature);
	}
	// Create a heating profile by reversing the cooling profile
	std::vector<double> heating_profile(temperature_profile.rbegin(), temperature_profile.rend());

	// Append the heating profile to the temperature profile
	temperature_profile.insert(temperature_profile.end(), heating_profile.begin(), heating_profile.end());
	return temperature_profile;
}

std::vector<double> record_occupancy(const std::map<int, std::vector<int>>& state) {
	std::vector<double> occupancy_vector;

	for (const auto& staple_state : state) {
		const std::vector<int>& domains = staple_state.second;

		std::map<int, int> counts;
		int total_domains = domains.size();

		for (int domain_state : domains) {
			if (domain_state != 0) {
				counts[domain_state]++;
			}
		}

		int max_count = 0;
		for (const auto& count : counts) {
			if (count.second > max_count) {
				max_count = count.second;
			}
		}

		double occupancy = static_cast<double>(max_count) / static_cast<double>(total_domains);
		occupancy_vector.push_back(occupancy);
	}

	return occupancy_vector;
}

std::vector<std::vector<int>> record_state(const std::map<int, std::vector<int>>& state) {
	std::vector<std::vector<int>> result;

	for (const auto& entry : state) {
		result.push_back(entry.second);
	}

	return result;
}

std::string state_to_string(const std::map<int, std::vector<int>>& state) {
	std::stringstream ss;

	for (const auto& staple : state) {
		for (const auto& domain : staple.second) {
			ss << domain;
		}
	}

	return ss.str();
}

std::vector<std::vector<std::vector<std::vector<int>>>> run_gillespie(std::vector<double>& temperature_profile, std::map<int, std::vector<std::vector<int>>>& stapPaths, std::map<int, std::vector<int>>& con,
	std::map<int, std::vector<std::vector<int>>>& vertex,
	std::map<int, std::vector<int>>& state,
	std::map<int, std::vector<std::vector<std::vector<int>>>>& neighbor,
	std::map<int, std::vector<std::vector<int>>>& sequence,
	std::vector<std::vector<double>>& adjacency_matrix,
	double time_interval,
	int number_trajectories, bool seq) {

	std::vector<std::vector<std::vector<std::vector<int>>>> recorded_states(number_trajectories);
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis(0.0, 1.0);

	for (int traj = 0; traj < number_trajectories; traj++) {
		std::map<int, std::vector<int>> state_copy = state;
		std::map<int, std::vector<std::vector<int>>> vertex_copy = vertex;
		std::map<int, std::vector<std::vector<std::vector<int>>>> neighbor_copy = neighbor;
		std::vector<std::vector<double>> adjacency_matrix_copy = adjacency_matrix;

		std::vector<std::vector<double>> trajectory_occupancies;
		std::vector<std::vector<std::vector<int>>> trajectory_states;

		// Add the initial state record for the first temperature
		if (!temperature_profile.empty()) {
			trajectory_states.push_back(record_state(state_copy));
		}

		// Start the loop from the first temperature
		for (auto it = std::next(temperature_profile.begin()); it != temperature_profile.end(); ++it) {
			double T = *it;
			double average_binding_energy = calculate_average_binding_energy(T);
			std::map<int, std::vector<double>> binding_energy = calculate_binding_energies(sequence, stapPaths, T, average_binding_energy, seq);

			double current_time = 0;

			while (current_time < time_interval) {
				std::pair<std::vector<Event>, double> events_propensity = handle_events(stapPaths, con, vertex_copy, state_copy, neighbor_copy, binding_energy, adjacency_matrix_copy, T, average_binding_energy);
				std::vector<Event>& events = events_propensity.first;
				double propensity = events_propensity.second;
				double delta_t = std::log(1.0 / dis(gen)) / propensity;
				current_time += delta_t;

				if (current_time < time_interval) {
					double r2 = dis(gen);
					double cumulative_rate = 0.0;

					for (const Event& event : events) {
						cumulative_rate += event.rate;
						if (r2 < cumulative_rate / propensity) {
							apply_changes(event, state_copy, adjacency_matrix_copy);
							break;
						}
					}
				}
			}

			trajectory_states.push_back(record_state(state_copy));
		}

		recorded_states[traj] = trajectory_states;

	}

	return recorded_states;
}

std::pair<std::vector<double>, std::vector<std::vector<std::vector<std::vector<int>>>>> run_gillespie_constant_T(
	double temperature,
	std::map<int, std::vector<std::vector<int>>>& stapPaths,
	std::map<int, std::vector<int>>& con,
	std::map<int, std::vector<std::vector<int>>>& vertex,
	std::map<int, std::vector<int>>& state,
	std::map<int, std::vector<std::vector<std::vector<int>>>>& neighbor,
	std::map<int, std::vector<std::vector<int>>>& sequence,
	std::vector<std::vector<double>>& adjacency_matrix,
	double time_interval,
	int number_trajectories, bool seq, double record_interval) {

	std::vector<std::vector<std::vector<std::vector<int>>>> recorded_states(number_trajectories);
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis(0.0, 1.0);
	auto start_time = std::chrono::steady_clock::now();
	std::chrono::seconds max_duration(71 * 60 * 60); // 71 hours
	std::vector<double> recorded_times_corrected;

	for (int traj = 0; traj < number_trajectories; traj++) {
		std::map<int, std::vector<int>> state_copy = state;
		std::map<int, std::vector<std::vector<int>>> vertex_copy = vertex;
		std::map<int, std::vector<std::vector<std::vector<int>>>> neighbor_copy = neighbor;
		std::vector<std::vector<double>> adjacency_matrix_copy = adjacency_matrix;

		std::vector<std::vector<std::vector<int>>> trajectory_states;

		// Record initial state at time 0
		trajectory_states.push_back(record_state(state_copy));
		recorded_times_corrected.push_back(0.0);

		double average_binding_energy = calculate_average_binding_energy(temperature);
		std::map<int, std::vector<double>> binding_energy = calculate_binding_energies(sequence, stapPaths, temperature, average_binding_energy, seq);

		double current_time = 0;
		double next_record_time = record_interval;
		while (current_time < time_interval) {
			//std::cout << current_time << std::endl;
			auto now = std::chrono::steady_clock::now();
			if (std::chrono::duration_cast<std::chrono::seconds>(now - start_time) > max_duration) {
				//std::cout << "Maximum duration reached. Stopping simulation..." << std::endl;
				recorded_states[traj] = trajectory_states;
				return { recorded_times_corrected, recorded_states };
			}
			std::pair<std::vector<Event>, double> events_propensity = handle_events(stapPaths, con, vertex_copy, state_copy, neighbor_copy, binding_energy, adjacency_matrix_copy, temperature, average_binding_energy);
			std::vector<Event>& events = events_propensity.first;
			double propensity = events_propensity.second;
			double delta_t = std::log(1.0 / dis(gen)) / propensity;
			//std::cout << current_time + delta_t << std::endl;
			if (current_time + delta_t < next_record_time) {
				// Select and perform event, then move time forward
				double r2 = dis(gen);
				double cumulative_rate = 0.0;
				Event selected_event;
				for (const Event& event : events) {
					cumulative_rate += event.rate;
					if (r2 < cumulative_rate / propensity) {
						selected_event = event;
						break;
					}
				}
				apply_changes(selected_event, state_copy, adjacency_matrix_copy);
				current_time += delta_t;
			}

			else if (current_time + delta_t > next_record_time) {
				// Record the current state and time, update current_time and next_record_time
				trajectory_states.push_back(record_state(state_copy));
				recorded_times_corrected.push_back(next_record_time);
				current_time = next_record_time;
				//std::cout << "We change the time to " << current_time << std::endl;
				//std::cout << state_to_string(state_copy) << std::endl;
				next_record_time += record_interval;
			}
		}
		recorded_states[traj] = trajectory_states;
	}

	return { recorded_times_corrected, recorded_states };
}



void print_recorded_states(const std::vector<std::vector<std::vector<double>>>& recorded_states) {
	for (size_t traj = 0; traj < recorded_states.size(); ++traj) {
		std::cout << "Trajectory " << traj + 1 << ":" << std::endl;

		const std::vector<std::vector<double>>& trajectory_states = recorded_states[traj];

		for (size_t temp_index = 0; temp_index < trajectory_states.size(); ++temp_index) {
			std::cout << "Temperature " << temp_index + 1 << " states:" << std::endl;

			const std::vector<double>& states = trajectory_states[temp_index];

			for (size_t state_index = 0; state_index < states.size(); ++state_index) {
				std::cout << "Staple " << state_index + 1 << ": " << states[state_index] << std::endl;
			}
		}

		std::cout << std::endl;
	}
}

std::string vector_to_string(const std::vector<int>& vec) {
	std::stringstream ss;
	for (size_t i = 0; i < vec.size(); ++i) {
		if (i != 0) {
			ss << "-";
		}
		ss << vec[i];
	}
	return ss.str();
}

void save_trajectory_states_to_csv(const std::vector<std::vector<std::vector<std::vector<int>>>>& recorded_state, const std::vector<double>& temperature_profile, const std::string& file_name) {
	int num_staples = recorded_state[0][0].size();
	int num_temperatures = temperature_profile.size();
	int num_trajectories = recorded_state.size();

	std::ofstream output_file(file_name);

	// Write the header row
	output_file << "Temperature";
	for (int traj_idx = 0; traj_idx < num_trajectories; traj_idx++) {
		for (int staple_idx = 0; staple_idx < num_staples; staple_idx++) {
			output_file << ",Trajectory_" << traj_idx << "_Staple_" << staple_idx;
		}
	}
	output_file << std::endl;

	// Write the data rows
	for (int temp_idx = 0; temp_idx < num_temperatures; temp_idx++) {
		output_file << temperature_profile[temp_idx];
		for (int traj_idx = 0; traj_idx < num_trajectories; traj_idx++) {
			for (int staple_idx = 0; staple_idx < num_staples; staple_idx++) {
				output_file << "," << vector_to_string(recorded_state[traj_idx][temp_idx][staple_idx]);
			}
		}
		output_file << std::endl;
	}

	output_file.close();
}
void save_trajectory_states_to_csv_T(const std::vector<std::vector<std::vector<std::vector<int>>>>& recorded_states, const std::vector<double>& recorded_times, const std::string& file_name) {
	if (recorded_states.empty() || recorded_states[0].empty() || recorded_times.empty()) {
		std::cerr << "Error: Empty data cannot be saved to CSV." << std::endl;
		return;
	}

	int num_staples = recorded_states[0][0].size();
	int num_times = recorded_times.size();
	int num_trajectories = recorded_states.size();

	std::ofstream output_file(file_name);

	// Write the header row
	output_file << "Time";
	for (int traj_idx = 0; traj_idx < num_trajectories; traj_idx++) {
		for (int staple_idx = 0; staple_idx < num_staples; staple_idx++) {
			output_file << ",Trajectory_" << traj_idx << "_Staple_" << staple_idx;
		}
	}
	output_file << std::endl;

	// Write the data rows
	for (int time_idx = 0; time_idx < num_times; time_idx++) {
		output_file << recorded_times[time_idx];
		for (int traj_idx = 0; traj_idx < num_trajectories; traj_idx++) {
			// Check if time_idx is within bounds for recorded_states[traj_idx]
			if (time_idx < recorded_states[traj_idx].size() && !recorded_states[traj_idx][time_idx].empty()) {
				for (int staple_idx = 0; staple_idx < num_staples; staple_idx++) {
					// Check if staple_idx is within bounds for recorded_states[traj_idx][time_idx]
					if (staple_idx < recorded_states[traj_idx][time_idx].size()) {
						output_file << "," << vector_to_string(recorded_states[traj_idx][time_idx][staple_idx]);
					}
				}
			}
		}
		output_file << std::endl;
	}

	output_file.close();
}

void print_recorded_states(const std::vector<double>& recorded_times, const std::vector<std::vector<std::vector<std::vector<int>>>>& recorded_states) {
	int num_trajectories = recorded_states.size();

	for (int traj_idx = 0; traj_idx < num_trajectories; traj_idx++) {
		std::cout << "Trajectory " << traj_idx << ":\n";

		for (size_t time_idx = 0; time_idx < recorded_times.size(); time_idx++) {
			std::cout << "Time " << std::fixed << std::setprecision(2) << recorded_times[time_idx] << ":\n";

			const auto& state = recorded_states[traj_idx][time_idx];
			for (size_t staple_idx = 0; staple_idx < state.size(); staple_idx++) {
				std::cout << "Staple " << staple_idx << ": ";

				for (int value : state[staple_idx]) {
					std::cout << value << " ";
				}

				std::cout << std::endl;
			}

			std::cout << std::endl;
		}

		std::cout << std::endl;
	}
}
void print_binding_energies(std::map<int, std::vector<std::vector<int>>> sequence, std::map<int, std::vector<std::vector<int>>> stapPaths, double temperature, double average_binding_energy, bool seq) {
	std::map<int, std::vector<double>> binding_energies = calculate_binding_energies(sequence, stapPaths, temperature, average_binding_energy, seq);

	for (const auto& staple_pair : binding_energies) {
		int staple_id = staple_pair.first;
		const auto& domain_binding_energies = staple_pair.second;

		double total_staple_energy = 0;

		std::cout << "Staple ID: " << staple_id << std::endl;
		std::cout << "Domain binding energies:" << std::endl;
		for (size_t i = 0; i < domain_binding_energies.size(); i++) {
			double energy = domain_binding_energies[i];
			total_staple_energy += energy;
			std::cout << "  Domain " << i << ": " << energy << std::endl;
		}

		double mean_energy = total_staple_energy / domain_binding_energies.size();
		double variance = 0;

		for (size_t i = 0; i < domain_binding_energies.size(); i++) {
			double diff = domain_binding_energies[i] - mean_energy;
			variance += diff * diff;
		}

		variance /= domain_binding_energies.size();
		double energy_stddev = std::sqrt(variance);

		std::cout << "Total staple binding energy: " << total_staple_energy << std::endl;
		std::cout << "Standard deviation of binding energy: " << energy_stddev << std::endl << std::endl;
	}
}
void save_binding_energies(std::map<int, std::vector<std::vector<int>>> sequence,
	std::map<int, std::vector<std::vector<int>>> stapPaths,
	double temperature,
	double average_binding_energy,
	bool seq) {

	std::map<int, std::vector<double>> binding_energies = calculate_binding_energies(sequence, stapPaths, temperature, average_binding_energy, seq);

	std::ofstream csv_file;
	csv_file.open("staple_binding_energies.csv");

	csv_file << "Staple ID, Total Staple Binding Energy\n";

	for (const auto& staple_pair : binding_energies) {
		int staple_id = staple_pair.first;
		const auto& domain_binding_energies = staple_pair.second;

		double total_staple_energy = 0;

		for (size_t i = 0; i < domain_binding_energies.size(); i++) {
			double energy = domain_binding_energies[i];
			total_staple_energy += energy;
		}

		csv_file << staple_id << "," << total_staple_energy << "\n";
	}

	csv_file.close();
}


std::vector<double> calculate_dGStaple(std::map<int, std::vector<int>>& con,
	std::map<int, std::vector<std::vector<int>>>& vertex,
	std::map<int, std::vector<std::vector<int>>>& stapPaths,
	std::vector<std::vector<double>>& adjacency_matrix,double T, double gamma, double w_c) {

	std::vector<double> dGShapeStaple;

	for (auto& staple : con) { // Loop through each staple in the connectivity map
		int staple_id = staple.first;
		std::vector<int>& connectivity = staple.second;

		double dGShape = 0.0;

		for (size_t i = 0; i < connectivity.size() - 1; i++) { // Loop through the connectivity of the current staple
			int v0 =0;
			int v1=0;

			switch (connectivity[i]) {
			case 0:
				v0 = vertex[staple_id][i][0];
				v1 = vertex[staple_id][i + 1][0];
				break;
			case 1:
				v0 = vertex[staple_id][i][1];
				v1 = vertex[staple_id][i + 1][1];
				break;
			case 2:
				v0 = vertex[staple_id][i][1];
				v1 = vertex[staple_id][i + 1][0];
				break;
			default:
				v0 = vertex[staple_id][i][0];
				v1 = vertex[staple_id][i + 1][1];
				break;
			}

			double shortest_path = dijkstra_shortest_path(adjacency_matrix, v0, v1);
			double E = shortest_path + w_c;
			double dGshape_staple = -R * T * gamma * std::log(C / E);

			dGShape += dGshape_staple;
		}

		dGShapeStaple.push_back(dGShape);
	}

	return dGShapeStaple;
}
void print_total_shape_energy(const std::vector<double>& dGShapeStaple) {
	for (size_t i = 0; i < dGShapeStaple.size(); ++i) {
		std::cout << "Total shape energy for staple " << i << ": " << dGShapeStaple[i] << std::endl;
	}
}


void calculate_and_save_dGStaple(std::map<int, std::vector<int>>& con,
	std::map<int, std::vector<std::vector<int>>>& vertex,
	std::map<int, std::vector<std::vector<int>>>& stapPaths,
	std::vector<std::vector<double>>& adjacency_matrix, double T, double gamma, double w_c) {

	std::ofstream csv_file;
	csv_file.open("staple_shape_energies.csv");

	csv_file << "Staple ID, Total Staple Shape Energy\n";

	for (auto& staple : con) { // Loop through each staple in the connectivity map
		int staple_id = staple.first;
		std::vector<int>& connectivity = staple.second;

		double dGShape = 0.0;

		for (size_t i = 0; i < connectivity.size() - 1; i++) { // Loop through the connectivity of the current staple
			int v0 = 0;
			int v1 = 0;

			switch (connectivity[i]) {
			case 0:
				v0 = vertex[staple_id][i][0];
				v1 = vertex[staple_id][i + 1][0];
				break;
			case 1:
				v0 = vertex[staple_id][i][1];
				v1 = vertex[staple_id][i + 1][1];
				break;
			case 2:
				v0 = vertex[staple_id][i][1];
				v1 = vertex[staple_id][i + 1][0];
				break;
			default:
				v0 = vertex[staple_id][i][0];
				v1 = vertex[staple_id][i + 1][1];
				break;
			}

			double shortest_path = dijkstra_shortest_path(adjacency_matrix, v0, v1);
			double E = shortest_path + w_c;
			double dGshape_staple = -R * T * gamma * std::log(C / E);

			dGShape += dGshape_staple;
		}

		csv_file << staple_id << "," << dGShape << "\n";
	}

	csv_file.close();
}

int main() {

std::map<int, std::vector<std::vector<int>>> stapPaths = {
  {0, {{2, 0, 15}, {1, 0, 15}, {0,0,15} }},
  {1, {{0, 16, 31}, {1, 16, 31}, {2,16,31} }},
  {2, {{2, 32, 47}, {1, 32, 47}, {0, 32,47} }},
  {3, {{0, 48, 63}, {1, 48, 63}, {2, 48,63} }},
  {4, {{2, 64, 79}, {1, 64, 79}, {0, 64,79} }},
  {5, {{0, 80, 95}, {1, 80, 95}, {2,80,95} }},
  {6, {{2, 96, 111}, {1, 96, 111}, {0,96,111} }},
  {7, {{0, 112, 127}, {1, 112, 127}, {2,112,127} }},
  {8, {{2, 128, 143}, {1, 128, 143}, {0,128,143} }} };

std::map<int, std::vector<int>> con = {
	  {0, {0, 1, 0}},
	  {1, {0, 1, 0}},
	  {2, {0, 1, 0}},
	  {3, {0, 1, 0}},
	  {4, {0, 1, 0}},
	  {5, {0, 1, 0}},
	  {6, {0, 1, 0}},
	  {7, {0, 1, 0}},
	  {8, {0, 1, 0}} };

std::map<int, std::vector<std::vector<int>>> vertex = {
	  {0, {{20, 21}, {10, 11},{0,1} }},
	  {1, {{1, 2}, {11, 12},{21,22}  }},
	  {2, {{22, 23}, {12, 13},{2,3} }},
	  {3, {{3, 4}, {13, 14}, {23, 24} }},
	  {4, {{24, 25}, {14, 15}, {4,5} }},
	  {5, {{5, 6}, {15, 16}, {25, 26}}},
	  {6, {{26, 27}, {16, 17}, {6, 7}}},
	  {7, {{7, 8}, {17, 18}, {27,28}}},
	  {8, {{28, 29}, {18, 19}, {8, 9}}} };

std::map<int, std::vector<std::vector<int>>> sequence = {
	 {0, {{2, 1, 4, 4, 1, 4, 2, 4, 4, 1, 4, 3, 2, 4, 2, 2}, {3, 1, 4, 3, 4, 2, 4, 3, 1, 2, 1, 2, 3, 2, 4, 3},{4, 1, 2, 4, 4, 1, 3, 2, 2, 1, 2, 4, 3, 1, 4, 1} }},
	 {1, {{1, 2, 4, 3, 3, 2, 3, 1, 3, 3, 2, 3, 4, 4, 3, 2}, {4, 2, 2, 2, 2, 1, 3, 4, 2, 2, 1, 4, 1, 1, 1, 3},{4, 4, 3, 3, 3, 4, 3, 2, 4, 4, 3, 1, 1, 3, 4, 4} }},
	 {2, {{4, 1, 4, 4, 3, 4, 1, 2, 1, 4, 4, 3, 3, 3, 3, 1}, {2, 3, 2, 3, 4, 2, 1, 3, 2, 2, 1, 2, 3, 3, 3, 2},{2, 4, 4, 4, 3, 3, 3, 1, 3, 1, 3, 2, 1, 4, 1, 3} }},
	 {3, {{1, 4, 1, 4, 4, 2, 2, 4, 1, 2, 3, 1, 4, 4, 2, 3}, {1, 2, 2, 2, 4, 3, 3, 3, 1, 2, 1, 4, 2, 1, 4, 4},{1, 2, 1, 3, 4, 2, 1, 3, 1, 2, 2, 3, 3, 3, 2, 2} }},
		 {4, {{4, 2, 1, 3, 3, 3, 3, 2, 3, 4, 4, 3, 3, 4, 3, 1}, {4, 3, 4, 3, 4, 3, 3, 4, 2, 3, 4, 3, 1, 4, 2, 1},{1, 4, 3, 2, 4, 1, 4, 4, 3, 1, 4, 1, 1, 3, 3, 1} }},
		 {5, {{2, 3, 2, 3, 2, 4, 1, 1, 4, 4, 2, 1, 1, 4, 1, 1}, {3, 1, 2, 4, 1, 2, 1, 2, 2, 2, 4, 2, 2, 3, 4, 3},{4, 3, 3, 1, 2, 3, 3, 3, 1, 2, 4, 2, 3, 3, 1, 3} }},
		 {6, {{3, 3, 2, 1, 2, 2, 3, 1, 3, 3, 1, 4, 1, 1, 2, 2}, {3, 3, 4, 2, 2, 3, 2, 2, 1, 1, 4, 4, 1, 2, 2, 4},{4, 3, 3, 4, 1, 3, 3, 1, 4, 4, 3, 3, 4, 4, 1, 3} }},
		 {7, {{3, 1, 4, 3, 3, 3, 2, 1, 3, 3, 4, 4, 4, 3, 3, 1}, {3, 3, 2, 4, 2, 2, 3, 1, 3, 4, 1, 3, 2, 3, 3, 1},{1, 1, 4, 3, 4, 1, 3, 1, 2, 1, 4, 4, 3, 4, 1, 2} }},
		 {8, {{4, 4, 4, 2, 3, 2, 2, 2, 4, 1, 3, 4, 1, 2, 3, 4}, {3, 2, 1, 1, 3, 2, 2, 2, 1, 2, 2, 2, 3, 4, 1, 4},{1, 1, 3, 2, 2, 1, 2, 3, 3, 1, 3, 1, 1, 3, 1, 1} }} };

std::map<int, std::vector<int>> state = {
	  {0, {0, 0, 0}},
	  {1, {0, 0, 0}},
	  {2, {0, 0, 0}},
	  {3, {0, 0, 0}},
	  {4, {0, 0, 0}},
	  {5, {0, 0, 0}},
	  {6, {0, 0, 0}},
	  {7, {0, 0, 0}},
	  {8, {0, 0, 0}} };

std::map<int, std::vector<std::vector<std::vector<int>>>> neighbor = {
	  {0, {{{1, 2}}, {{1, 1}},{{1,0}} }},
	  {1, {{{0, 2}, {2, 2}}, {{0, 1}, {2, 1}},{{0,0}, {2,0}} }},
	  {2, {{{1, 2}, {3, 2}}, {{1, 1}, {3, 1}},{{1,0}, {3,0}} }},
	  {3, {{{2, 2}, {4, 2}}, {{2, 1}, {4, 1}},{{2,0}, {4,0}} }},
	  {4, {{{3, 2}, {5, 2}}, {{3, 1}, {5, 1}},{{3,0}, {5,0}} }},
	  {5, {{{4, 2}, {6, 2}}, {{4, 1}, {6, 1}},{{4,0}, {6,0}} }},
	  {6, {{{5, 2}, {7, 2}}, {{5, 1}, {7, 1}},{{5,0}, {7,0}} }},
	  {7, {{{6, 2}, {8, 2}}, {{6, 1}, {8, 1}},{{6,0}, {8,0}} }},
	  {8, {{{7, 2}}, {{7, 1}},{{7,0}} }} };

std::map<int, std::vector<int>> scafcross = {
	  {0, {0}},
	  {1, {143}} };

std::map<std::vector<int>, double> dangling_ends = {};
std::vector<std::vector<double>> adjacency_matrix = create_adjacency_matrix(stapPaths, vertex, scafcross, dangling_ends);
save_binding_energies(sequence, stapPaths, 323.15, 0.0, seq);
calculate_and_save_dGStaple(con, vertex, stapPaths, adjacency_matrix, 323.15, gamma_val, w_c);
//std::pair<std::vector<double>, std::vector<std::vector<std::vector<std::vector<int>>>>> recorded_info = run_gillespie_constant_T(328.15, stapPaths, con, vertex, state, neighbor, sequence, adjacency_matrix, 360, n_traj, seq, 0.05);
//save_trajectory_states_to_csv_T(recorded_info.second, recorded_info.first, "Constant_55.csv");
std::vector<double> temp_ramp = generate_temp_ramp(T_max, T_min, dTdt);
std::vector<std::vector<std::vector<std::vector<int>>>> recorded_info = run_gillespie(temp_ramp, stapPaths, con, vertex, state, neighbor, sequence, adjacency_matrix, dt, n_traj, seq);
save_trajectory_states_to_csv(recorded_info, temp_ramp, "5.csv");

	return 0;
}
