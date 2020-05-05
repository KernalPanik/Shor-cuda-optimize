#include <cstdlib>
#include "quantumRegister.h"
#include "rand.h"
#include "util.h"
#include <iostream>

#include <chrono>

#include "qft.h"

// Register class

Register::Register(unsigned int num_qubits) {
	/*
	By default, initializes vaccum (|000...>) to amplitude
	one and the rest zero. This can be adjust with logic
	gates, of course, or explicitly with the function
	set_nonzero_states.
	*/
	
	this->num_qubits = num_qubits;
	states[string(num_qubits, '0')] = amp(1, 0); // total prob 1;
}

vec_states Register::all_states(unsigned int n) {
	/* 
	returns vector of all the possible n qubit states
	IN ORDER, where we choose the basis to be in increasing
	order in terms of their binary representation.
	i.e. 0..00 < 0..01 < 0..10 < 0..11 < 1..00 < 1..01 < 1..11
	*/
	vec_states v;
	if (n == 1) {v.push_back("0"); v.push_back("1");}
	else if(n > 1) {
		for (string s : all_states(n - 1)) {
			v.push_back(s + "0");
			v.push_back(s + "1");
		}
	}
	return v;
}

bool Register::check_state(string state) {
	// See if state is in state map.
	return states.count(state);
}

void Register::set_nonzero_states(state_map &s) {
	amp total = 0.0;
	for (state_map::iterator i = s.begin(); i != s.end(); ++i) { 
		states[i->first] = i->second; total += pow(abs(i->second), 2); 
	}
	if (total == 0.0)
		printf("Bad input: must have at least one state with nonzero amplitude\n");
	else if (total == 1.0) 
		states = s;
	else { // normalize input
		states = s;
		for (state_map::iterator i = states.begin(); i != states.end(); ++i) 
			states[i->first] = i->second / sqrt(total);
	};
}

amp Register::amplitude(string state) {
	// Note: not actually physically measureable.
	if(check_state(state)) return states[state];
	else return 0;
}

double Register::probability(string state) {
	return pow(abs(amplitude(state)), 2);
}

string Register::measure() {
	/*
	Measure the system, and collapse it into a state.
	Update the states map to have 1 probability of being
	in this state, and return the state.
	*/
	// Will always return something because there is no way for the total
	// probability to not equal 1.
	double r = get_rand();
	double total(0.0); string s;
	for (state_map::iterator i = states.begin(); i != states.end(); ++i) {
		s = i->first;
		total += probability(s);
		if (r <= total) { // collapse
			// get rid of all other states which just take up memory with amp zero.
			states = state_map();
			states[s] = 1.0;
			return s;
		}
	}

	// WE SHOULD NEVER MAKE IT OUT HERE. If you use my functions, I make sure
    // that the total probability is always 1, so we will never make it out here.
    // But, if you adjusted the states manually with the operator[] funcionality
    // that I included and you did something wrong, then we may make it out here.
    // This is why I recommend you never use the reg[] functionality, but if you do,
    // make sure you get the physics correct so that the total probability is always
    // one and we never make it out here. In the event that we do, just pretend like
    // we measured |000...>.
    s = string(num_qubits, '0');
    states = state_map();
    states[s] = 1.0;
    return s;
    
}

char Register::measure(unsigned int qubit) {
	/*
	Measure a qubit, and collapse it.
	Update the states map to have 1 probability of having a particular
	value for the qubit, and returns the value.
	*/
	// Will always return something because there is no way for the total
	// probability to not equal 1.

	double zero_prob = 0; string s;
	for (state_map::iterator i = states.begin(); i != states.end(); ++i) {
		s = i->first;
		if(s[qubit] == '0') zero_prob += probability(s);
	}

	state_map m; char v; amp p;
	if (get_rand() < zero_prob) {v = '0'; p = sqrt(zero_prob);}
	else { v = '1'; p = sqrt(1 - zero_prob); }
	// Resete state map to only states with the measured qubit.
	for (state_map::iterator i = states.begin(); i != states.end(); ++i) {
		s = i->first;
		if (s[qubit] == v) m[s] = states[s] / p;
	}
	states = m;
	return v;
}

void Register::print_states() {
	// Show all nonzero states.
	cout << *this;
}

std::ostream &operator<<(std::ostream &os, Register &reg) {
	// Show all nonzero states
	for (state_map::iterator i = reg.states.begin(); i != reg.states.end(); ++i)
		os << "|" << i->first << ">: amp-> "
		   << (i->second).real() << " + " << (i->second).imag() << "i"
		   << ", prob-> " << reg.probability(i->first) << "\n";
	return os;
}

state_map Register::copy_map(state_map &s) {
	state_map m;
	for (state_map::iterator i = s.begin(); i != s.end(); ++i) { m[i->first] = i->second; }
	return m;
}

vec_states Register::nonzero_states() {
	/*
	Returns all the keys in the states map; so only
	states with nonzero amplitude.

	This function is primarily here to be used with the function
	below (operator[]). I HEAVILY RECOMMEND NOT USING THIS OR
	THAT FUNCTION. See below for why.
	*/
	vec_states v;
	for (state_map::iterator i = states.begin(); i != states.end(); ++i)
		v.push_back(i->first);
	return v;
}

amp & Register::operator[](string state) {
	/*
	Gives public access to the private states map. We return a reference,
	so that one can change the states dictionary from outside the class.
	THERE ARE NO CHECKS ON THIS FUNCTION. Anything you change, YOU must
	make sure that it acts as a unitary transformation, otherwise probabilites
	will fail, and then you can have issues with the measure function.

	I HEAVILY RECOMMEND NOT USING THIS FUNCTIONALITY. I make the states
	map private because you should really only be using quantum gates
	(unitary transformations) to affect the register, and should not be
	directly accessing the states. With that being said, I include this
	function as a last option.

	Use: 
		Register reg(4);
		cout << reg["0000"] << endl;
		reg["0000"] = 1 / sqrt(2);
		reg["1000"] = 1 / sqrt(2);
		cout << reg << endl;
	*/
	return states[state];
}


// Gates

void Register::apply_cufft(vec_int qubits)
{
	string state, s;
	unsigned int r, j;
	state_map old = copy_map(states);
	amp c;
	vec_states temp_states = all_states(qubits.size());
	std::vector<string> statenames;
	std::vector<amp> statevals;
	for (state_map::iterator i = old.begin(); i != old.end(); ++i) 
	{
		statenames.push_back(i->first);
		statevals.push_back(i->second);
	}

	auto newmap = CalculateQFT(statenames, statevals);
}

void Register::apply_gate(Unitary u, vec_int qubits) {
	/*
	Applys the unitary matrix u to the given qubits in the system.
	To get rid of unneccessary memory, if we come across a state with
	amplitude zero, remove it from the states map.

	Example:
		if vec_int qubits = [0, 2], then expect u.dimension to be 4.
		We apply the matrix u to the zeroth and second qubits in the
		system. The matrix is represented in the basis 
			{|00>, |01>, |10>, |11>}
		where the first number applies to the zeroth qubit and the
		second number applies to the second qubit.

	Example:
		if vec_int qubits = [2, 4, 0], then we expect u.dimension to be 8.
		We apply the matrix u to the second, fourth, and zeroth qubits in the
		system. The matrix is represented in the basis
			{|000>, |001>, |010>, |011>, |100>, |101>, |110>, |111>}
		where the first number applies to the second qubit, the second
		number applies the fourth qubit, and the third number applies
		to the zeroth qubit.

	All integers in vec_qubits should be different! I do not perform a check,
	but instead assume that the user uses the gates correctly.
	*/
	if (u.dimension != (unsigned int)(1 << qubits.size())) { // 1 << qubits.size is pow(2, qubits.size())
		printf("Unitary matrix dimension is not correct to be applied to the inputs qubits\n");
		return; 
	}

	string state, s;
	unsigned int r, j; 
	state_map old = copy_map(states); 
	amp c;
	vec_states temp_states = all_states(qubits.size());
	for (state_map::iterator i = old.begin(); i != old.end(); ++i) {
		state = i->first; 
		s = "";
		for (unsigned int q : qubits)
			s += state[q];

		r = binary_to_base10(s); // Find which number basis element s corresponds to.
		states[state] -= (1.0 - u[r][r]) * old[state];
		// if (states[state] == 0.0) states.erase(state); // Get rid of it.
		if (probability(state) < 1e-16) 
			states.erase(state); // zero beyond machine precision.
		j = 0;
		for(string k : temp_states) {
			if (j != r) {
				s = state;
				for (unsigned int l = 0; l < k.size(); l++) 
					s[qubits[l]] = k[l];
				c = u[j][r] * old[state];
				if (check_state(s)) {
					states[s] += c;
					// if (states[s] == 0.0) states.erase(s);
					if (probability(s) < 1e-16)
					 states.erase(s); // zero beyond machine precision.
				} else if(c != 0.0) 
					states[s] = c;
			}
			j++;
		}
	}
}

// Common gates listed below. Any gate you want to use that is not
// listed below can be implemented by just creating the unitary operator
// corresponding to the gate and calling the apply_gate function.

void Register::Hadamard(unsigned int qubit) {
	/*
	zero -> n-1 qubit indexing.
	Hadamard operator on single qubit is
		H = ((|0> + |1>) <0| + (|0> - |1>) <1|) / sqrt(2)
	*/

	vec_int v; 
    v.push_back(qubit);
	apply_gate(Unitary::Hadamard(), v);
}

void Register::PhaseShift(unsigned int qubit, double theta) {
	/*
	zero -> n-1 qubit indexing.
	Phase shift by theta is
		P = |1><0| + exp(i theta) |0><1|
	*/
	vec_int v; v.push_back(qubit);
	apply_gate(Unitary::PhaseShift(theta), v);
}

void Register::PiOverEight(unsigned int qubit) {
	// zero->n - 1 qubit indexing.
	// PhaseShift(qubit, pi / 4.0);
	vec_int v; v.push_back(qubit);
	apply_gate(Unitary::PiOverEight(), v);
}

void Register::PauliX(unsigned int qubit) {
	/*
	zero->n - 1 qubit indexing.
	Pauli-X gates, also the NOT gate, for a single qubit is
		PX = |1><0| + |0><1|
	*/
	vec_int v; v.push_back(qubit);
	apply_gate(Unitary::PauliX(), v);
}

void Register::PauliY(unsigned int qubit) {
	/*
	zero->n - 1 qubit indexing.
	Pauli-Y gate for a single qubit is
		PY = i(|1><0| - |0><1|)
	*/
	vec_int v; v.push_back(qubit);
	apply_gate(Unitary::PauliY(), v);
}

void Register::PauliZ(unsigned int qubit) {
	/*
	zero->n - 1 qubit indexing.
	Pauli-Z gate for a single qubit is
		PZ = |1><0| - |0><1|
	*/
	// PhaseShift(qubit, pi);
	vec_int v; v.push_back(qubit);
	apply_gate(Unitary::PauliZ(), v);
}

void Register::ControlledNot(unsigned int control_qubit, unsigned int target_qubit) {
	/*
	zero -> num_qubits-1 qubit indexing.
	ControlledNot gate is just the NOT gate (PauliX) on the target 
	qubit if the controlled qubit is 1. Otherwise, do nothing.
	*/
	vec_int v; v.push_back(control_qubit); v.push_back(target_qubit);
	apply_gate(Unitary::ControlledNot(), v);
}

void Register::Toffoli(unsigned int control_qubit1, unsigned int control_qubit2, unsigned int target_qubit) {
	/*
	zero -> num_qubits-1 qubit indexing.
	Toffoli gate, also known as the Controlled-Controlled-Not gate 
	is just the NOT gate (PauliX) on the target qubit if both the
	controlled qubits are 1. Otherwise, do nothing.
	*/
	vec_int v; v.push_back(control_qubit1); v.push_back(control_qubit2); v.push_back(target_qubit);
	apply_gate(Unitary::Toffoli(), v);
}

void Register::ControlledPhaseShift(unsigned int control_qubit, unsigned int target_qubit, double theta) {
	/*
	zero -> num_qubits-1 qubit indexing.
	Just the phase shift gate on the target qubit if the first qubit is 1.
	*/
	vec_int v; v.push_back(control_qubit); v.push_back(target_qubit);
	apply_gate(Unitary::ControlledPhaseShift(theta), v);
}

void Register::Swap(unsigned int qubit1, unsigned int qubit2) {
	/*
	zero -> num_qubits-1 qubit indexing.
	Swap qubit1 and qubit2.
	*/

	// vec_int v; v.push_back(qubit1); v.push_back(qubit2);
	// apply_gate(Unitary::Swap(), v);
	ControlledNot(qubit1, qubit2);
	ControlledNot(qubit2, qubit1);
	ControlledNot(qubit1, qubit2);
}

void Register::Ising(unsigned int qubit1, unsigned int qubit2, double theta) {
	vec_int v; v.push_back(qubit1); v.push_back(qubit2);
	apply_gate(Unitary::Ising(theta), v);
}


// Sort of a gate
void Register::apply_function(function<string(string)> f) {
	/*
	Unitary transformation that sends 
		|x> to |f(x)>

	This one is not technically a gate in the code, but in real life would be
	a unitary transformation that could be implemented with a quantum circuit.

	Note that in order for this to be unitary, f(x) must map each x to a unique
	f(x). ie. f(x) cannot equal f(y) unless x = y. This is checked in this function.
	*/
	
	state_map old = copy_map(states); string state; states = state_map();
	for (state_map::iterator i = old.begin(); i != old.end(); ++i) {
		state = f(i->first);

		if (check_state(state)) {
			printf("The provided function does not constitute a unitary transformation.\n");
			states = old; return; // reset.
		}

		states[state] = i->second;
	}
	
}
