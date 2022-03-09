function seq = ReorderNoReps(seq)

  % Find unique values and counts:
  N = numel(seq);
  [vals, ~, index] = unique(seq(:), 'stable');
  counts = accumarray(index, 1);
  [maxCount, maxIndex] = max(counts);

  % Check the maximum number of occurrences:
  if (2*maxCount-1 > N)
    error('Can''t produce sequence without repeats!');
  end

  % Fill cell array column-wise with permuted and replicated set of values:
  C = cell(maxCount, ceil(N/maxCount));
  if (3*maxCount-1 > N)
    permIndex = [maxIndex(1) ...
                 setdiff(randperm(numel(vals)), maxIndex(1), 'stable')];
  else
    permIndex = randperm(numel(vals));
  end
  C(1:N) = num2cell(repelem(vals(permIndex), counts(permIndex)));

  % Transpose cell array and extract non-empty entries:
  C = C.';
  seq = reshape([C{~cellfun('isempty', C)}], size(seq));

end