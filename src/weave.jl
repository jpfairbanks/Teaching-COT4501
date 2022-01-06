using Weave

files = ["epimodel.jmd"]
for file in files
    weave(file)
    tangle(file)
end
