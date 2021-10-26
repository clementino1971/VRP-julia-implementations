using JuMP,CPLEX

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

    c_pos = [Vector{Int64}(undef,2) for _ in 1:N]
    q = [0 for _ in 1:N]
    for i in 8:8+N-1
        str = split(lines[i]," ")
        c_pos[parse(Int64,str[2])] = [parse(Int64,str[3]),parse(Int64,str[4])]
    end 
    
    for i in 8+N+1:8+N+1+N-1
        str = split(lines[i]," ")
        q[parse(Int64,str[1])] = parse(Int64,str[2]) 
    end 


    return instance(c_pos,q,1,k,Q,N)
end

function run_cvrp(G)
    CVRP = Model(with_optimizer(CPLEX.Optimizer))
    @variable(CVRP,x[1:G.N,1:G.N],Bin)
    @variable(CVRP,u[1:G.N],lower_bound = 0, upper_bound = G.Q)

    #=States that in a route, each customer vertex is connected 
    to two other vertices
    =# 
    for i in 2:G.N
        @constraint(CVRP, sum(x[i,1:i-1]) + sum(x[i,i+1:G.N]) == 1)
        @constraint(CVRP, sum(x[1:i-1,i]) + sum(x[i+1:G.N,i]) == 1)
    end

    #=
        Ensures that at most k routes are constructed
    =#
    @constraint(CVRP, sum(x[1,1:G.N]) == G.k)


    # Miller-Tucker-Zemlin inequalities
    for i=2:G.N
        for j=2:G.N
            if i != j
                @constraint(CVRP,u[i] - u[j] + G.Q*x[i,j] <= G.Q-G.q[j])
            end
        end
    end

    for i in 2:G.N
        @constraint(CVRP,G.q[i] <= u[i] <= G.Q)
    end

    @objective(CVRP,Min, sum(x[i,j]*distance(G.customers[i],G.customers[j]) for i=1:G.N,j=1:G.N))

    println(CVRP)

    optimize!(CVRP)
    @show solution_summary(CVRP, verbose=true)
end

function main(args)
    if length(args) != 1
        prinln("Wrong args number! The call shoud be cvrp-mtz.jl <instance-filename>") 
    else
        G = read_instance(args[1])
        run_cvrp(G)
    end
end

main(ARGS)