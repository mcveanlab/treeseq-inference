import tempfile
import subprocess
import logging
import string

import msprime, pyslim

eidos_cmd = string.Template("""
function (void)save_treeseq(string prefix, string string_freq, integer ogens) {
	catn("Saving in generation " + sim.generation);
	sim.treeSeqOutput(prefix+string_freq+(ogens ? format("+%i", ogens) else "")+".decap");
}

function (void)check_freq_save(integer freq_index) {
	if (sim.countOfMutationsOfType(m1)) {
		//there is a segregating site
		if (sim.mutationFrequencies(NULL, NULL) > output_at_freq_float[freq_index]) {
			//we have hit the first freq
			save_treeseq(treefile_prefix, output_at_freq[freq_index], 0);
			// don't need this event any more
			sim.deregisterScriptBlock(freq_index);
			// but we might need to check the next freq in line
			if ((freq_index+1)<length(output_at_freq_float)) {
				register_freq_check(freq_index+1);
			} else {
				// No more frequencies to check
				if (length(output_after_fixation_gens)==0) sim.simulationFinished();
			}
		}
	}
}

function (void)register_freq_check(integer f_index) {
	sim.registerLateEvent(f_index, format("{check_freq_save(%i);}", f_index), $equilibration_gens, NULL);
}

initialize() {
    $set_random_seed_cmd
	initializeTreeSeq();
	initializeMutationRate(0);
	initializeMutationType("m1", $dominance_coefficient, "f", $selection_coefficient);
	initializeGenomicElementType("g1", m1, 1.0);
	initializeGenomicElement(g1, 0, $length-1);
	initializeRecombinationRate($recombination_rate);
	output_at_frequency = c($freq_strings); // strings, for nice filename creation
	output_generations_after_fixation = c($output_gens);
	// Save these variables into constants so we can access them in other functions
	freq_float = asFloat(output_at_frequency);
	defineConstant("output_at_freq",output_at_frequency[order(freq_float)]);
	defineConstant("output_at_freq_float", freq_float[order(freq_float)]);
	defineConstant("output_after_fixation_gens", output_generations_after_fixation);
	defineConstant("treefile_prefix", "$treefile_prefix");
}

1 late() {
    sim.addSubpop("p0", $popsize);
	register_freq_check(0);
}

// Allow {equilibration_gens} generations for equilibrating before 
$equilibration_gens: late() {
	if (sim.substitutions.size()) {
		catn("Fixed in generation " + sim.generation);
		// Now fixed: stop checking for fixation / reintroduction
		sim.deregisterScriptBlock(self);
		for (post_gen in output_after_fixation_gens) {
			event = format("save_treeseq(treefile_prefix, '1.0', %i);", post_gen);
			if (post_gen == max(output_after_fixation_gens)) 
				event = event+"sim.simulationFinished();";
			if (post_gen == 0) {
				executeLambda(event); // Should output NOW
			} else {
				output_gen = post_gen + sim.generation; // Output in further generations
				sim.registerLateEvent(NULL, "{"+event+"}", output_gen, output_gen);
			}
		}
	} else if (sim.countOfMutationsOfType(m1) == 0) { //no mutations, must introduce
		target = sample(sim.subpopulations.genomes, 1);
		target.addNewDrawnMutation(m1, $mutant_position);
		catn("Introduced in generation " + sim.generation);
	}
}

$max_generations {
	stop("Got to end without fixing");
}
""")

def comma_separated_list(arr):
    '''
    Return string version of list without braces, e.g. for [1,2,3] return "1,2,3" and
    for ['a','b','c'] return "'a','b','c'". This should work for both arrays, sets, and
    tuples (since e.g. sets are formatted as {1,2,3}
    '''
    return "{}".format(arr)[1:-1]


def recapitate_amd_mutate(filename, mu, rho, Ne, seed):
    ts = pyslim.load(filename) #no simplify
    ts = ts.recapitate(recombination_rate=rho, Ne=Ne, random_seed=seed)
    return msprime.mutate(ts, mu, random_seed=seed, keep=True)


def simulate_sweep(popsize, chrom_length, recomb_rate, mut_rate, 
    selection_coef, dominance_coef, nsamples, output_at_freqs, 
    mutations_after_simulation = True, equilibration_gens=100,
    max_generations=1e9, treefile_prefix="sweepfile", seed=None, slimname="slim"):
    """
    Carry out a simulation of a selective sweep, and save msprime-format files at frequencies
    specified by output_at_freqs, which is a list of (frequency, post_generation) tuples
    Note that some of these files may have fixed variants (i.e. mutations above the root node)
    
    nsamples is number of *diploid* samples
    """
    freq_to_output = set()
    gens_post_fixation_to_output = set()
    for o in output_at_freqs:
        is_tuple = isinstance(o, tuple)
        freq = o[0] if is_tuple else o
        post_gens = int(o[1]) if is_tuple else 0
        assert isinstance(freq, str)
        if float(freq) == 1.0:
            #this is specified as something at fixation
            gens_post_fixation_to_output.add(post_gens)
        else:
            assert post_gens==0 # Ban output some generations after an intermediate freq
            freq_to_output.add(freq)

    eidos_seed_cmd = "" if seed is None else "setSeed({});".format(int(seed))
    with tempfile.NamedTemporaryFile(mode="w+t") as eidos_file:
        #use a temporary file for the SLiM commands until https://github.com/MesserLab/SLiM/issues/24 is addressed 
        cmd = eidos_cmd.substitute(
            set_random_seed_cmd = eidos_seed_cmd,
            treefile_prefix = treefile_prefix,
            dominance_coefficient = dominance_coef,
            selection_coefficient = selection_coef,
            mutant_position = chrom_length//2,
            popsize = popsize,
            length = chrom_length,
            recombination_rate = recomb_rate,
            freq_strings = comma_separated_list(freq_to_output),
            output_gens  = comma_separated_list(gens_post_fixation_to_output),
            max_generations = int(max_generations),
            equilibration_gens = equilibration_gens
        )
        print(cmd, file=eidos_file, flush=True)
        subprocess.call([slimname,eidos_file.name], stdout=subprocess.DEVNULL)
    ret_val = {}
    for o in output_at_freqs:
        is_tuple = isinstance(o, tuple)
        freq = o[0] if is_tuple else o        
        fn = treefile_prefix + freq
        if is_tuple and len(o)>1 and o[1]:
            fn += "+%i" % o[1]
        ts = recapitate_amd_mutate(fn + ".decap", mut_rate, recomb_rate, popsize, seed)
        fn += ".trees"
        ts.dump(fn)
        ret_val[o] = fn
    return ret_val

if __name__ == '__main__':
    simulate_sweep(
        5000, 100000, 1e-8, 1e-8, 0.1, 0.5, 16,  
        output_at_freqs=['0.2', '0.5', '0.8', '1.0', ('1.0', 200), ('1.0', 100)])