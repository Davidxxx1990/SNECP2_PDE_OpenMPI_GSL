function plot_uxt5_t8(file_t5, file_t8)
	ux_t5 = load(file_t5);
	ux_t8 = load(file_t8);

	plot(ux_t5(:,1), ux_t5(:,2)); hold on;
	plot(ux_t8(:,1), ux_t8(:,2));
	
	grid on;
	xlabel('t');
	ylabel('exitation');
	xlabel('space');

	xlim([0 0.5]);
	
	legend('u(x,5)', 'u(x,8)');
	
end
