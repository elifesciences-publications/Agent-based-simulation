# Agent-based-simulation
Agent based simulation of constant and active disassembly processes to simulate length regulation in Chlamydomonas.

Each program takes parameters which the user can change-: time over which the process is simulated in minutes(n) Diffusion constant(D), total tubulin in micron(T), gamma in the model(gammma), value of gamma*k*M(gammakM), constant disassembly constant term(da), velovity of ballistic motion(v), steady state length(Lss), time at which flagella is severed(in minutes), time step(dt), (replenishment time for motors and tubulins(t_m and t_t respectively);only  in TM_shared_depoly)

Files include-:

TM_shared_depoly-: Agent based simulation of active disassembly process with the option of including feedback and severing

