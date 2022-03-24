function greedy_navigate(source, destination, graph, distm; verbose = true)
    success = true
    cur = source
    prev = -1
    route = [cur]
    i = 1
    while cur != destination && i <= size(distm, 1)
        i += 1
        possible_next = setdiff(graph.weights[:, cur].nzind, cur)
        if length(possible_next) == 0
            if verbose 
                print("Cannot navigate from $source to $destination: ")
                println("$cur is an absorbing node.")
            end
            success = false
            break
        end
        next = possible_next[findmin(pn -> distm[pn, destination], possible_next)[2]]
        if next == prev
            if verbose 
                print("Cannot navigate from $source to $destination: ")
                println("$cur and $next form a greedy loop.")
            end
            success = false
            break
        end
        append!(route, next)
        prev, cur = cur, next
    end
    return route, success
end

function navigability(digraph, dm, N = 1000)
    n = size(digraph.weights, 1)
    success_rate = 0.0
    mean_length = 0.0
    for i in 1:N
        s, t = sample(1:n, 2, replace = false)
        route, success = greedy_navigate(s, t, digraph, dm, verbose = false)
        lr = length(route)
        success_val = Int(success)
        success_rate = success_rate + (success_val - success_rate)/i
        mean_length = mean_length + (lr - mean_length)/i
    end
    mean_length, success_rate
end