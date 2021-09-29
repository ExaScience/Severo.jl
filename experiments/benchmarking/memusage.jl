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

struct RUsage
    ru_utime_sec::Clong         #  user CPU time used
    ru_utime_usec::Clong        #  user CPU time used
    ru_stime_sec::Clong         #  system CPU time used
    ru_stime_usec::Clong        #  system CPU time used
    ru_maxrss::Clong            #  maximum resident set size
    ru_ixrss::Clong             #  integral shared memory size
    ru_idrss::Clong             #  integral unshared data size
    ru_isrss::Clong             #  integral unshared stack size
    ru_minflt::Clong            #  page reclaims (soft page faults)
    ru_majflt::Clong            #  page faults (hard page faults)
    ru_nswap::Clong             #  swaps
    ru_inblock::Clong           #  block input operations
    ru_oublock::Clong           #  block output operations
    ru_msgsnd::Clong            #  IPC messages sent
    ru_msgrcv::Clong            #  IPC messages received
    ru_nsignals::Clong          #  signals received
    ru_nvcsw::Clong             #  voluntary context switches
    ru_nivcsw::Clong            #  involuntary context switches
end

function get_vmsize()
    ru = Ref{RUsage}()
    ccall(:getrusage, Cint, (Cint, Ptr{Cvoid}), 0, ru)
    return ru[].ru_maxrss
end

