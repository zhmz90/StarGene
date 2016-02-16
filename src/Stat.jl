module Stat

const cosmic_fl = "../data/CosmicCompleteExport.tsv"

const data_dir = "../data"

function heredgenes()
    s = Set(readcsv("data/58gene.csv",ASCIIString))
    push!(s,"BRCA1")
    push!(s,"BRCA2")
    s
end


const panel_fl = "/home/guo/haplox/Bone/data/panel/Cancer_gene.list"

function genespanel()
    genespanel(panel_fl)
end

@doc """ Get sorted genes list from a panel file .list
""" ->
function genespanel(fl)
    genes = Set{ASCIIString}()
    open(fl) do file
        header = readline(file)
        while !eof(file)
            line = readline(file)
            gene = strip(line,'\n')
            if !in(gene,genes)
                push!(genes,gene)
            end
        end
    end
    filter!(gene->gene!="KAT6A" && gene!="MEK1",genes)
    push!(genes,"MYST3")
    push!(genes,"MAP2K1")
    sort(collect(genes))
end

function check_heredgene()
    hdgenes = heredgenes()
    cgenes  = cosmicgenes()
    sygns = Set{ASCIIString}()
    for gn in hdgenes
        if in(gn,cgenes)
            continue
        end
        push!(sygns,gn)
    end
    sygns
end

function check_pnlhdg()
    data = collect(setdiff(Set(genespanel()), heredgenes()))
    writecsv("panel_genes_not_hereditary.csv",sort(data))
end
function cosmicgenes()
    genes = Set{ASCIIString}()
    open(cosmic_fl) do file
        header = readline(file)
        while !eof(file)
            line = readline(file)
            row = split(line,"\t")
            gene = row[1]
            if contains(gene,"_")
                gene = split(gene,"_")[1]
            end
            gene = convert(ASCIIString,gene)
            push!(genes,gene)
        end
    end
    genes
end
function check_pnlgenes()
    pnlgenes = genespanel()
    cosgenes = cosmicgenes()
    for gene in pnlgenes
        if !in(gene,cosgenes)
            @show gene
        end
    end
end

function cancernames()
    cancers = Set{ASCIIString}()
#    cancer_gene = Dict{ASCIIString,Dict{ASCIIString,Arr}}()
    open(cosmic_fl) do file
        header = readline(file)
#        @show split(header,"\t")[17]
#        @show header
        while !eof(file)
            line = readline(file)
            row = split(line,"\t")
 #           gene = convert(ASCIIString, row[1])
            cancer = convert(ASCIIString, row[8])
 #           sample_name = convert(ASCIIString,row[5])
 #           mut = length(row[17]) != 0
 #           @show gene,cancer,sample_name,mut
 #           break
            push!(cancers,cancer)
        end
    end
    collect(cancers)
end

function genesforcancer(cancer)
    gene_samples = Dict{ASCIIString,Dict{ASCIIString,Bool}}()
    samples = Set{ASCIIString}()
    #hdgenes = heredgenes()
    pnlgenes = genespanel()
    open(cosmic_fl) do file
        header = readline(file)
        while !eof(file)
            line = readline(file)
            row = split(line,"\t")
            cancer_name = convert(ASCIIString, row[8])
            if cancer_name != cancer
                continue
            end
            sample_name = convert(ASCIIString, row[5])
            if !in(sample_name, samples)
                push!(samples, sample_name)
            end
            gene = convert(ASCIIString, row[1])
            if contains(gene,"_")
                gene = convert(ASCIIString, split(gene,"_")[1])
            end
            if !in(gene, pnlgenes)
                continue
            end
            mut = length(row[17]) != 0
            if !in(gene, keys(gene_samples))
                gene_samples[gene] = Dict{ASCIIString,Bool}(sample_name=>mut)
            else
                if in(sample_name, keys(gene_samples[gene]))
                    gene_samples[gene][sample_name] |= mut
                else
                    gene_samples[gene][sample_name] = mut
                end
            end
        end
    end
    genes = collect(keys(gene_samples)) #first column genename
    samplemuts = map(collect, collect(values(gene_samples))) # [[sp1=>0,sp2=>1],[],[]]
    numsamples = map(length, samplemuts) # third column tol num samples
    nummuts = map(samplemuts) do spmut # second column num mutations
        length(filter(x->x[2]==1,spmut))
    end

    ratios = nummuts ./ numsamples # fourth column of mutation rates

    #=
    samplemuts = map(length,collect(values(gene_samples)))
    data = sort(collect(Dict(zip(genes,counts))),by=x->x[2],rev=true)
    gene = map(x->x[1],data)
    len  = length(samples)
    count = map(x->x[2],data)
    ratio = map(x->x[2]/float(len),data)
    lens  = repmat([len],size(data,1),1)
    =#
    
    header = ["Gene name","The number of muated samples",
              "The total number of samples in $cancer","mutation rates"]
    
    data = hcat(genes, nummuts, numsamples, ratios)
    data = sortrows(data,by=x->(x[3],x[2]),rev=true)
    data = data[(data[:,4] .>= 0.1) & (data[:,3] .>= 10),:]
    view = vcat(header', data)
    writecsv(joinpath(data_dir, "$(cancer).csv"), view)
end

function genesampleofeachcancer()
    cancers = cancernames()
    pmap(genesforcancer, cancers)
end


end
