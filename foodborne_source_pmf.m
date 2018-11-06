% --- INPUTS ---
% flows = The Markov transition matrix for the directed, layered network G
% stage_ends(stage) = last node in layer/stage.
% prior_pmf = row vect of prior probability such that prior_pmf(s) = probability of source node s if contam_observed_nodes is empty
% contam_reports = reported locations of illness. first row is the node_ID contamination observed at; second row is the corresponding time
% Frac = Fraction of contaminated nodes a source node needs to reach to be considered feasible

% --- OUTPUTS ---
% pmf = a row vector of the probability of each feasible source being the true source

function [pmf] = foodborne_source_pmf(flows, stage_ends, prior_pmf, contam_reports, Frac)

pmf = zeros([1, stage_ends(1)]); % pmf will eventually be row vector of a probability for each feasible source
flows = sparse(flows); % for runtime

% identify feasible sources as those which reach a fraction >= Frac of all contaminated nodes
feasible_sources = 1:stage_ends(1); % initialize to list of all farm nodes
reaches = zeros(size(feasible_sources)); % contains num of contaminated nodes reached by farm of same index
for c_node = unique(contam_reports(1,:)) % compute all potential sources
    last_parents = c_node;
    while ~isempty(parents(last_parents, flows))
        last_parents = parents(last_parents, flows);
    end % end while
    reaches(last_parents) = reaches(last_parents) + 1; % update reach counts
end % end for
feasible_sources = find(reaches >= Frac*length(unique(contam_reports(1,:)))); % feasible sources IDs are those which reach a fraction-greater-than-P of contaminated nodes

% compute pmf using Markov chain matrix method, combine with prior pmf, and normalize to sum to 1
pmf_from_network = pmf_over_sources(feasible_sources, contam_reports, flows, stage_ends);
pmf = pmf_from_network.*prior_pmf; % elementwise multiply
pmf = pmf / sum(pmf); % normalize

end % end function


% Given a list of feasible sources, the Markov transition matrix, and stage_ends, 
% calculate the probability distribution of feasible sources to be the true source,
% return the pmf over sources.
% Formula is A = (I_R - P_Q)^{-1}*P_R
function pmf = pmf_over_sources(feasible_sources, contam_reports, flows, stage_ends)
    pmf = zeros([1, stage_ends(1)]);
    
    % Extract submatrices from Markov transition matrix 
    num_stages = length(stage_ends); % for convenience
    I_R = eye(stage_ends(num_stages-1)); %submatrix of size |V_R| X |V_R| representing absorption at an absorbing node (i.e. a node in last stage)
    P_Q = flows(1:stage_ends(end-1), 1:stage_ends(end-1)); %submatrix of size |V_Q| X |V_Q| concerning transitions between transient nodes (i.e. nodes in all stages but the last)
    P_R = flows(1:stage_ends(end-1), stage_ends(end-1)+1:stage_ends(end)); %submatrix of size |V_Q| X |V_R| concerning transitions from transient nodes into absorbing nodes

    % Calculate absorbing probability matrix
    A = inv(I_R - P_Q) * P_R;

    % Calculate probability of source s generating the whole collection of contamination reports
    for s = feasible_sources
        diff_traj_likelihood = 1; % multiply all path probabilities together onto this variable       
        for node_ID = contam_reports(1, :)
           path_likelihood = A(s, node_ID-stage_ends(end-1)); % likelihood that the current observation started at s, summed over all possible paths
           diff_traj_likelihood = diff_traj_likelihood * path_likelihood; % probability of source s generating the whole collection of reports
        end % end for       
        pmf(s) = diff_traj_likelihood;
    end 
end

% node_set = array of node IDs that we are querying for the overall set of ancestors of. Assumed to all be in same stage.
% flows = the adjacency matrix
% prents = an array containing all the IDs of parents of node_set
function prents = parents(node_set, flows)
    prents = [];
    for j = 1:length(node_set)
        prents = union(prents, find(flows(:, node_set(j)))); % find each node that leads into nodes in j
        [r, c] = size(prents);
        if r > c % make sure it comes out as a row vector
            prents = prents.';
        end % end if
    end % end for
end % end ancestors
