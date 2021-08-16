function d = get_design(test, reader, id)
  counts2D = crosstab(test, reader);
  counts3D = crosstab(test, reader, id);
  
  d.balanced2D = all(counts2D, 'all');
  d.balanced3D = all(counts3D == 1, 'all');

  f = @(x, y) all(sum(crosstab(y, x) > 0, 1) == 1, 'all');
  d.id_in_reader = d.balanced2D && f(id, reader);
  d.id_in_test = d.balanced2D && f(id, test);
  d.reader_in_test = f(reader, test);
end
