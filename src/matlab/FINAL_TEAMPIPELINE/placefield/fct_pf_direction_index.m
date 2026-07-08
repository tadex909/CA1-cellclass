function [direction_index] = fct_pf_direction_index(firing_rate_way1, firing_rate_way2)


direction_index = abs(nansum(firing_rate_way1 - firing_rate_way2(end:-1:1)) / nansum(firing_rate_way1 + firing_rate_way2(end:-1:1)));