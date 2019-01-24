function si = calculate_similarity_index(trl,rdm,stim)

% get stim value
sv = stim(trl);

% change stimulus value if retrieval
if trl > size(rdm,1)/2
    sv = sv - 4;
end

%  find trls of same value
match_idx = find(stim==sv);
match_idx(match_idx==trl) = []; % <- remove identical trial

% find mismatching trials at encoding
mismatch_idx = find(stim(stim<=4)~=sv);

% calculate similarity index
si = -(mean(rdm(trl,match_idx)) - mean(rdm(trl,mismatch_idx))); %#ok<FNDSB>
