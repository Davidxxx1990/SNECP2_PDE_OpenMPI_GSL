function plot_ut3L4_L2(file_3L4, file_L2)
	ut_3L4 = load(file_3L4);
	ut_L2 = load(file_L2);

	plot(ut_3L4(:,1), ut_3L4(:,2)); hold on;
	plot(ut_L2(:,1), ut_L2(:,2));
	
	grid on;
	xlabel('t');
	ylabel('exitation');
	xlabel('time');

	xlim([0 10]);
	
	legend('u(x,5)', 'u(x,8)');
	
end
