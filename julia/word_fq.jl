filename = ARGS[1]
customlt(a,b) = (b.second < a.second) ? true : b.second == a.second ? a.first < b.first : false

function wordfreq(filename::String)
  file = open(filename)
  freqtab = Dict{String, Int64}()
  while !eof(file)
    for word in split(readline(file))
      index = Base.ht_keyindex(freqtab, word)
      if index > 0
        freqtab.vals[index] += 1
      else
        freqtab[word] = 1
      end
    end
  end
  v = collect(freqtab); sort!(v,lt=customlt)
  # v = collect(freqtab); sort!(v, alg=QuickSort)
  # v = collect(freqtab); sort!(v, alg=QuickSort, rev=true)
  # for t in v
  #   @printf("%s\t%d\n",t.first,t.second)
end

@time tab = wordfreq(ARGS[1])
# gc_enable(false)
