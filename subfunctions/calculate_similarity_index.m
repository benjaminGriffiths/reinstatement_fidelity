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
X = [rdm(trl,match_idx) rdm(trl,mismatch_idx)];
Y = [zeros(1,numel(match_idx))-1 zeros(1,numel(mismatch_idx))];
si = atanh(corr(X',Y'));
