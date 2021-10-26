using JuMP,CPLEX,LightGraphs,GLPK,LinearAlgebra

struct instance
    customers
    q
    depot
    k
    Q
    N
end

function distance(u,v)
    
    x_sq = (v[1] - u[1])^2
    y_sq = (v[2] - u[2])^2
    #if data.round
    return floor(sqrt(x_sq + y_sq) + 0.5)
    #end
    #return sqrt(x_sq + y_sq)
end
 
# This function read instances on format of CVRPLIB available at: http://vrp.atd-lab.inf.puc-rio.br/index.php/en/
# and transform this data on a instance struct above
function read_instance(file_name)
    println("filename: ",file_name)
    f = open(file_name)
    lines = readlines(f)
    for i in 1:3
        println(lines[i])
    end


    aux = last(split(lines[1],"-"))
    k = parse(Int64,aux[2:length(aux)])
    N = parse(Int64, split(lines[4]," ")[3])
    Q = parse(Int64, split(lines[6]," ")[3])

    c_pos = [Vector{Float64}(undef,2) for _ in 1:N]
    q = [0 for _ in 1:N]
    for i in 8:8+N-1
        str = split(lines[i]," ")
        #println(str)
        c_pos[parse(Int64,str[2])] = [parse(Float64,str[3]),parse(Float64,str[4])]
    end 
    
    for i in 8+N+1:8+N+1+N-1
        str = split(lines[i]," ")
        #println(str)
        q[parse(Int64,str[1])] = parse(Int64,str[2]) 
    end 


    return instance(c_pos,q,1,k,Q,N)
end

function delta(inst, S, x)
    L = []

    for i in S
        for j in 1:inst.N
            if !(j in S)
                i_min = min(i,j)
                i_max = max(i,j)
                if !(x[i_min,i_max] in L)
                    append!(L, [x[i_min,i_max]])
                end
            end
        end
    end
    return L
end

function demand_s(inst, S)
    d = 0
    for c in S
        d = d + inst.q[c]
    end
    return d
end



function delta_cb(inst, S, x)
    sm = 0

    L = []

    for i in S
        for j in 1:inst.N
            if !(j in S)
                i_min = min(i,j)
                i_max = max(i,j)
                println(i_min," ",i_max, " " ,x[i_min][i_max])
                if (!([i_min,i_max] in L)) && (x[i_min][i_max] > 0.0)
                    sm += x[i_min][i_max]
                    append!(L, [[i_min,i_max]])
                end
            end
        end
    end

    println(L)
    return sm
end


function W(inst, x, curr_s)
    L = delta_cb(inst, curr_s, x)
    println(L)
    return L - 2*ceil(demand_s(inst, curr_s)/inst.Q)
end

function build_simple_graph(inst, x)
    g = SimpleGraph(inst.N)
    for i in 2:nv(g)
        for j in i+1:nv(g)
            if(value(x[i][j]) != 0)
                 add_edge!(g, i, j)
            end
        end
     end
     return g
end



function run_cvrp(G)
     
    CVRP = direct_model(CPLEX.Optimizer()) 
    set_silent(CVRP)

    MOI.set(CVRP, MOI.NumberOfThreads(), 1)

    x = @variable(CVRP,x[i=1:G.N,i+1:G.N],Int)
    for i in 1:G.N
        for j in i+1:G.N
            if(i == 1)
                c = @constraint(CVRP, 0 <= x[i,j])
                c = @constraint(CVRP, x[i,j] <= 2)
            else
                c = @constraint(CVRP, 0 <= x[i,j])
                c = @constraint(CVRP, x[i,j] <= 1)
            end;
        end;
    end;


    #=States that in a route, each customer vertex is connected 
    to two other vertices
    =# 
    for i in 2:G.N
        dt = delta(G, [i], x)
        if(length(dt) == 0)
            continue
        end

        c = @constraint(CVRP, sum(dt) == 2 )
    end

    #=
        Ensures that at most k routes are constructed
    =#
    
    c = @constraint(CVRP, sum(delta(G, [1], x)) == 2*G.k )

    @objective(CVRP,Min, sum(x[i,j]*distance(G.customers[i],G.customers[j]) for i=1:G.N,j=i+1:G.N))

    println(CVRP)

    cnt = 1
    cb_calls = Clong[]

    function ressource_constraints(cb_data::CPLEX.CallbackContext, context_id::Clong)
                
        # You can select where the callback is run
        if context_id != CPX_CALLBACKCONTEXT_CANDIDATE
            return
        end

        #println("\n #### SOLUTION ", cnt, "####")
        cnt += 1
        
        CPLEX.load_callback_variable_primal(cb_data, context_id)
        nx = [[0.0 for _ in 1:G.N] for _ in 1:G.N]
        #pre = [callback_value(cb_data,x_k ) for x_k in x]
        
        k = 1
        for i in 2:G.N
            for j in i+1:G.N
                nx[i][j] = callback_value(cb_data,x[i,j])
            end
        end
        
        g = build_simple_graph(G, nx)
        comp = connected_components(g)

        for j in 2:G.N
            nx[1][j] = callback_value(cb_data,x[1,j])
        end
        
        for c in comp
            println("\nComponent: ", c)
            deleteat!(c, findall(x->x==1,c))
            println(nx)
            f = W(G, nx,c)
            println("f: ",f)
            if(f<0)
                dt = delta(G, c, x)
                if(length(dt) > 0)
                    con = @build_constraint(sum(dt) >= 2*ceil(demand_s(G, c)/G.Q))
 
                    println(con)

                    MOI.submit(CVRP, MOI.LazyConstraint(cb_data), con)
                end
            end
        end
    end

    MOI.set(CVRP, CPLEX.CallbackFunction(), ressource_constraints)
    
    optimize!(CVRP)
    @show solution_summary(CVRP, verbose=true)
end

function main(args)
    if length(args) != 1
        prinln("Wrong args number! The call shoud be cvrp-bc.jl <instance-filename>") 
    else
        G = read_instance(args[1])
        run_cvrp(G)
    end
end

main(ARGS)


