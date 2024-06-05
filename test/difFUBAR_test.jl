begin
    analysis_name = "nobackground/Ace2"
    seqnames,seqs = read_fasta("data/Ace2_no_background/Ace2_tiny_test.fasta");
    treestring, tags, tag_colors = import_colored_figtree_nexus_as_tagged_tree("data/Ace2_no_background/Ace2_no_background.nex")
    df,results = difFUBAR(seqnames, seqs, treestring, tags, tag_colors, analysis_name, exports = false, verbosity = 0)
    @assert size(df) == (19, 8)
end

begin #getpuresubclades
    function two_pure_clades!(tree::FelNode, tags::Vector{String})
        for node in getnodelist(tree.children[1])
            node.name *= first(tags)
        end
        for node in getnodelist(tree.children[2])
            node.name *= last(tags)
        end
    end
    
    function randomly_assign_tags!(tree::FelNode, tags::Vector{String})
        for node in getnodelist(tree)
            node.name *= rand(vcat(tags, ""))
        end
    end    

    #Two pure clades
    tree = sim_tree(n=20)
    tags = ["{G1}", "{G2}"]
    two_pure_clades!(tree, tags)
    pure_subclades = CodonMolecularEvolution.getpuresubclades(tree, tags)
    @test length(pure_subclades) == 2 || any(map(isleafnode, tree.children))
    @test isroot(pure_subclades[1].parent) && isroot(pure_subclades[2].parent)

    #Random tags
    tree = sim_tree(n=200)
    tags = ["{G1}", "{G2}"]
    randomly_assign_tags!(tree, tags)
    pure_subclades = CodonMolecularEvolution.getpuresubclades(tree, tags)
    tag_ind(n) = CodonMolecularEvolution.model_ind(n.name, tags)
    @test all(pure_subclade -> begin
        nodelist_below = getnodelist(pure_subclade)[2:end]
        some_tag_ind = tag_ind(first(nodelist_below))
        tag_inds_below = map(tag_ind, nodelist_below)
        #Criteria for nodes below
        all(some_tag_ind .== tag_inds_below)
    end, pure_subclades)
end