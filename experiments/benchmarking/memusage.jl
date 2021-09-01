function get_memusage(pid=getpid())
    open("/proc/$pid/statm") do io
        line = readline(io)
        parts = split(line)
        vmrss = parse(Int64, parts[2])
        4096 * vmrss
    end
end

function memusage()
    l = ReentrantLock()
    peak_usage = 0

    t = @task begin
        while true
            try
                lock(l)
                peak_usage = max(peak_usage, get_memusage())
 #                ccall(:printf, Int, (Ptr{UInt8},Int32), "%d\n", peak_usage)
            finally
                unlock(l)
            end

            sleep(0.01)
        end
    end
    schedule(t)

    function query_memusage()
        rv = try
            lock(l)
            rv = peak_usage
            peak_usage = get_memusage()
            rv
        finally
            unlock(l)
        end
        println("mem $rv")
        rv
    end
end
