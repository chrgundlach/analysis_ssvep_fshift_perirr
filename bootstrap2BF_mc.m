function [BF,options] = bootstrap2BF_mc(dist1,dist2, plotflag, varargin)
% This function takes two bootstrapped distributions (>1000 draws needed)
% dist1 and dist2, which reflect an effect model (dist1) and a null model
% (dist2). It can compute a BF for a directional hypothesis (positive or negative)
% or for a meaningful vs. negligible difference (ROPE{region of practical equivalence} analysis). It computes 
% then the BF as posterior oddsover prior odds for the
% model entering the bootstrap in dist1. In many cases dist2 will be a
% permutation distribution representing the null model. 
% dist1 and dist2 are numeric vectors

% to avoid issues with scaling, do a joint z-normalization of both
% bootstrapped distributions together, thus maintaining any differences in
% mean and spread
%
% epsilon_orig is in units of the input (e.g. ms).
%
% input distributions a bootstrapped parameter estimate and it is advisable
% to have 10000 samples at least. to estimate the standard error an estimated proportion p is roughly sqrt(p(1−p)/N​).
%
% note: this is an adapted version of Andreas Keil's function
% bootstrap2BF_z.m
%
% 2025 - Norman Forschack
%
% --- Input Parsing ---
p = inputParser;
addRequired(p, 'dist1', @isnumeric);
addRequired(p, 'dist2', @isnumeric);
addRequired(p, 'plotflag', @(x) ismember(x, [0,1]));


% Define parameters for different test types
addParameter(p, 'test_type', 'directional_positive', ...
             @(x) ismember(x, {'directional_positive', 'directional_negative', 'bidirectional'}));
addParameter(p, 'epsilon_orig', NaN, ... % ROPE half-width in original units
             @(x) isnumeric(x) && isscalar(x)); % User must set this if test_type is 'rope'
parse(p, dist1, dist2, plotflag, varargin{:});
options = p.Results;

% --- Data Preparation ---
dist1_orig = dist1(:); % Keep original for ROPE calculation
dist2_orig = dist2(:); % Keep original for ROPE calculation

% Joint z-normalization
fulldist_orig = [dist1_orig; dist2_orig];
m_fd = mean(fulldist_orig,'omitnan');
s_fd = std(fulldist_orig,'omitnan');

% Handle cases with zero standard deviation to avoid division by zero
% MATLAB's zscore produces all zeros if std is zero. We'll mimic or rely on normalize.
if s_fd < sqrt(eps) % A very small number, effectively zero
    warning('Standard deviation of combined distributions is near zero. Z-scores may not be meaningful.');
    % normalize might produce NaNs or Zeros. If Zeros, directional test probabilites will be 0 or 1.
    % If NaNs, probabilities will be NaN.
    dist1z = (dist1_orig - m_fd); % Will be all zeros if s_fd = 0
    dist2z = (dist2_orig - m_fd); % Will be all zeros if s_fd = 0
    if s_fd ~=0 % if s_fd is small but non-zero, still divide.
       dist1z = dist1z ./ s_fd;
       dist2z = dist2z ./ s_fd;
    end
else
    dist1z = (dist1_orig - m_fd) ./ s_fd;
    dist2z = (dist2_orig - m_fd) ./ s_fd;
end

% Create versions of dist1z and dist2z with NaNs removed
dist1z_valid = dist1z(~isnan(dist1z));
dist2z_valid = dist2z(~isnan(dist2z));


% % Check for empty arrays after NaN removal to avoid division by zero
% if isempty(dist1z_valid)
%     posteriorsignedlikelyhood_effect_empirical = NaN; % Or 0, or handle as error
%     warning('effect probability cannot be estimated because all samples contained nan values!')
% else
%     posteriorsignedlikelyhood_effect_empirical = sum(dist1z_valid > 0) / length(dist1z_valid);
% end
% 
% if isempty(dist2z_valid)
%     signedlikelyhood_null_empirical = NaN; % Or 0, or handle as error
% else
%     signedlikelyhood_null_empirical = sum(dist2z_valid > 0) / length(dist2z_valid);
% end

% --- Calculate Probabilities Based on Test Type ---
prob_H1_effect_model = NaN; % Probability of the hypothesis of interest under the effect model
prob_H1_null_model = NaN;   % Probability of the hypothesis of interest under the null model

% Check for empty arrays after NaN removal to avoid division by zero (exit)
if isempty(dist1z_valid)
    error('effect probability cannot be estimated because all samples contained nan values!')
elseif isempty(dist2z_valid)
    error('effect probability under the null hypothesis cannot be estimated because all samples contained nan values!')
elseif isempty(dist1z_valid) && (dist2z_valid)
    error('the two input distributions only contain nan values. calculation of odds-ratios not possible!')
else
    switch options.test_type
        case 'directional_positive'
            % H1: difference > 0
            if ~isnan(options.epsilon_orig)
                prob_H1_effect_model = sum(dist1_orig > options.epsilon_orig) / length(dist1_orig);
                prob_H1_null_model   = sum(dist2_orig > options.epsilon_orig) / length(dist2_orig);
            else
                prob_H1_effect_model = sum(dist1z_valid > 0) / length(dist1z_valid);
                prob_H1_null_model   = sum(dist2z_valid > 0) / length(dist2z_valid);
            end
            
        case 'directional_negative'
            % H1: difference < 0
            if ~isnan(options.epsilon_orig)
                prob_H1_effect_model = sum(dist1_orig < options.epsilon_orig) / length(dist1_orig);
                prob_H1_null_model   = sum(dist2_orig < options.epsilon_orig) / length(dist2_orig);
            else
                prob_H1_effect_model = sum(dist1z_valid < 0) / length(dist1z_valid);
                prob_H1_null_model   = sum(dist2z_valid < 0) / length(dist2z_valid);
            end
            
        case 'bidirectional'
            if options.epsilon_orig <= 0
                error('For ROPE test, epsilon_orig (ROPE threshold in original units) must be a positive value specified by the user.');
            end
            % H1: difference is meaningful (i.e., |difference| > epsilon_z)
            if ~isnan(options.epsilon_orig)
                prob_H1_effect_model = sum(abs(dist1_orig) > options.epsilon_orig) / length(dist1_orig(~isnan(dist1_orig)));
                prob_H1_null_model   = sum(abs(dist2_orig) > options.epsilon_orig) / length(dist2_orig(~isnan(dist2_orig)));
            else
                prob_H1_effect_model = sum(abs(dist1z_valid) > 0) / length(dist1z_valid(~isnan(dist1z_valid)));
                prob_H1_null_model   = sum(abs(dist2z_valid) > 0) / length(dist2z_valid(~isnan(dist2z_valid)));
            end
            
        otherwise
            error('Unknown test_type specified.');
    end
end
% --- Calculate precision (standard error) of likelihoods 
LH_precision = [sum(isnan(dist1z)) sum(isnan(dist2z));... % absolute number of nans within each distribution
    sum(isnan(dist1z))/numel(dist1z) sum(isnan(dist2z))/numel(dist2z);... % relative number of nans
    sqrt((prob_H1_effect_model*(1-prob_H1_effect_model))/sum(~isnan(dist1z))) sqrt((prob_H1_null_model*(1-prob_H1_null_model))/sum(~isnan(dist2z)))]; % standard error of the estimated proportion (likelihood)

% --- Odds Calculation ---
% Odds = P(H1 is true) / P(H1 is false) = P(H1) / (1 - P(H1))

odds_posterior = prob_H1_effect_model / (1 - prob_H1_effect_model);
if prob_H1_effect_model == 1
    odds_posterior = Inf;
elseif prob_H1_effect_model == 0
    odds_posterior = 0;
end

odds_prior = prob_H1_null_model / (1 - prob_H1_null_model);
if prob_H1_null_model == 1
    odds_prior = Inf;
elseif prob_H1_null_model == 0
    odds_prior = 0;
end


% --- Bayes Factor Calculation ---
if odds_prior == 0 && odds_posterior == 0
    BF = 1; % Or NaN; 
    disp('both models agree H1 is impossible (a difference is unlikely under both hypotheses).')
elseif odds_prior == 0 && odds_posterior > 0
    BF = Inf; % 
    disp('Null model: H1 impossible; Effect model: H1 possible.')
elseif isinf(odds_prior) && isinf(odds_posterior)
    BF = 1; % Or NaN; 
    disp('both models agree H1 is certain (a difference is infinetely likely under both hypotheses).')
elseif isinf(odds_prior) && ~isinf(odds_posterior)
    BF = 0; % 
    disp('Null model: H1 certain; Effect model: H1 not certain.')
else
    BF = odds_posterior / odds_prior;
end

% --- Calculate and store additional statistics ---
% The CIs should be calculated on the original scale data
options.stats.ci_effect_95 = prctile(dist1_orig, [2.5, 97.5]);
options.stats.ci_null_95   = prctile(dist2_orig, [2.5, 97.5]);
options.stats.mean_effect  = mean(dist1_orig,'omitnan');
options.stats.median_effect= median(dist1_orig,'omitnan');
options.stats.sd_effect    = std(dist1_orig,'omitnan');
options.stats.mean_null    = mean(dist2_orig,'omitnan');
options.stats.median_null  = median(dist2_orig,'omitnan');
options.stats.sd_null      = std(dist2_orig,'omitnan');
options.stats.LH_precision = LH_precision; % store likelihood/fidelity of hypotheses
options.stats.odds_posterior_effect = odds_posterior;
options.stats.odds_posterior_null = odds_prior; 

% fit the bootstrapped distributions as normal distributions
% PD1 = fitdist(column(dist1z),'normal');
% PD2 = fitdist(column(dist2z),'normal');
% PD1 = fitdist(dist1z(:),'normal');
% PD2 = fitdist(dist2z(:),'normal');

% apply the normal distributions
% x_values = -5:.01:5;
% y1 = pdf(PD1,x_values);
% y2 = pdf(PD2,x_values);

% plot for control, then turn off plotflag
if plotflag == 1
%     figure;
% subplot(2,1,1), histogram(dist1z, 100, 'Normalization','pdf')
% title('empirical')
% hold on
% plot(x_values, y1, 'LineWidth',3)
% xline(0)
% xlim([-6 6])
% 
% subplot(2,1,2), histogram(dist2z, 100, 'Normalization','pdf')
% title('random')
% hold on
% plot(x_values, y2, 'LineWidth',3)
% xline(0)
% xlim([-6 6])

figure;
    x_hist_range = linspace(min([dist1z; dist2z]), max([dist1z; dist2z]), 101);

    subplot(2,1,1);
    histogram(dist1z, x_hist_range, 'Normalization','pdf');
    title(['Effect Model (dist1z) - Test: ', strrep(options.test_type, '_', ' ')]);
    hold on;
    if strcmp(options.test_type, 'bidirectional')
        if s_fd > sqrt(eps) && ~isnan(options.epsilon_orig) % Only draw ROPE on z-plot if s_fd is not effectively zero
            z_rope_upper_on_plot = (options.epsilon_orig - m_fd) / s_fd;
            z_rope_lower_on_plot = (-options.epsilon_orig - m_fd) / s_fd;
            
            actual_z_rope_lower_on_plot = min(z_rope_lower_on_plot, z_rope_upper_on_plot);
            actual_z_rope_upper_on_plot = max(z_rope_lower_on_plot, z_rope_upper_on_plot);

            % On each subplot (gca):
            y_lim = get(gca, 'YLim');
            fill([actual_z_rope_lower_on_plot, actual_z_rope_upper_on_plot, actual_z_rope_upper_on_plot, actual_z_rope_lower_on_plot],...
                [y_lim(1), y_lim(1), y_lim(2), y_lim(2)], ...
                'm', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'DisplayName', 'ROPE');
%             xline(actual_z_rope_lower_on_plot, 'm--'); xline(actual_z_rope_upper_on_plot, 'm--');
        else
            % Optionally, add text to plot indicating s_fd is ~0 and how to interpret original ROPE
%             text(pos_x, pos_y, sprintf('Original data std is ~0. Values ~%.2f. Original ROPE is +/-%.2f ms', m_fd, options.epsilon_orig));
        end
    elseif strcmp(options.test_type, 'directional_positive')
        if s_fd > sqrt(eps) && ~isnan(options.epsilon_orig) % Only draw ROPE on z-plot if s_fd is not effectively zero
            z_rope_upper_on_plot = (options.epsilon_orig - m_fd) / s_fd;
            z_rope_lower_on_plot = (0 - m_fd) / s_fd;

            actual_z_rope_lower_on_plot = min(z_rope_lower_on_plot, z_rope_upper_on_plot);
            actual_z_rope_upper_on_plot = max(z_rope_lower_on_plot, z_rope_upper_on_plot);

            y_lim = get(gca, 'YLim');
            fill([actual_z_rope_lower_on_plot, actual_z_rope_upper_on_plot, actual_z_rope_upper_on_plot, actual_z_rope_lower_on_plot],...
                [y_lim(1), y_lim(1), y_lim(2), y_lim(2)], ...
                'm', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'DisplayName', 'ROPE');
        else
            % Optionally, add text to plot indicating s_fd is ~0 and how to interpret original ROPE
%             text(pos_x, pos_y, sprintf('Original data std is ~0. Values ~%.2f. Original ROPE is +/-%.2f ms', m_fd, options.epsilon_orig));
        end

    elseif strcmp(options.test_type, 'directional_negative')

        if s_fd > sqrt(eps) && ~isnan(options.epsilon_orig) % Only draw ROPE on z-plot if s_fd is not effectively zero
            z_rope_upper_on_plot = (0 - m_fd) / s_fd;
            z_rope_lower_on_plot = (options.epsilon_orig - m_fd) / s_fd;

            actual_z_rope_lower_on_plot = min(z_rope_lower_on_plot, z_rope_upper_on_plot);
            actual_z_rope_upper_on_plot = max(z_rope_lower_on_plot, z_rope_upper_on_plot);

            y_lim = get(gca, 'YLim');
            fill([actual_z_rope_lower_on_plot, actual_z_rope_upper_on_plot, actual_z_rope_upper_on_plot, actual_z_rope_lower_on_plot],...
                [y_lim(1), y_lim(1), y_lim(2), y_lim(2)], ...
                'm', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'DisplayName', 'ROPE');
        else
            % Optionally, add text to plot indicating s_fd is ~0 and how to interpret original ROPE
%             text(pos_x, pos_y, sprintf('Original data std is ~0. Values ~%.2f. Original ROPE is +/-%.2f ms', m_fd, options.epsilon_orig));
        end

    else
        xline(0, 'k--');
    end
    xlim([min(x_hist_range) max(x_hist_range)]);
    xlabel('Normalized Latency Difference (z-score)');
    ylabel('Density');
    legend('show', 'Location', 'best'); % Show ROPE in legend if applicable
    hold off;

    subplot(2,1,2);
    histogram(dist2z, x_hist_range, 'Normalization','pdf');
    title(['Null Model (dist2z) - Test: ', strrep(options.test_type, '_', ' ')]);
    hold on;
    if strcmp(options.test_type, 'bidirectional')
        if s_fd > sqrt(eps)  && ~isnan(options.epsilon_orig) % Only draw ROPE on z-plot if s_fd is not effectively zero
            z_rope_upper_on_plot = (options.epsilon_orig - m_fd) / s_fd;
            z_rope_lower_on_plot = (-options.epsilon_orig - m_fd) / s_fd;
            
            actual_z_rope_lower_on_plot = min(z_rope_lower_on_plot, z_rope_upper_on_plot);
            actual_z_rope_upper_on_plot = max(z_rope_lower_on_plot, z_rope_upper_on_plot);

            % On each subplot (gca):
            y_lim = get(gca, 'YLim');
            fill([actual_z_rope_lower_on_plot, actual_z_rope_upper_on_plot, actual_z_rope_upper_on_plot, actual_z_rope_lower_on_plot],...
                [y_lim(1), y_lim(1), y_lim(2), y_lim(2)], ...
                'm', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'DisplayName', 'ROPE');
%             xline(actual_z_rope_lower_on_plot, 'm--'); xline(actual_z_rope_upper_on_plot, 'm--');
        else
            % Optionally, add text to plot indicating s_fd is ~0 and how to interpret original ROPE
%             text(pos_x, pos_y, sprintf('Original data std is ~0. Values ~%.2f. Original ROPE is +/-%.2f ms', m_fd, options.epsilon_orig));
        end
    elseif strcmp(options.test_type, 'directional_positive')
        if s_fd > sqrt(eps)  && ~isnan(options.epsilon_orig) % Only draw ROPE on z-plot if s_fd is not effectively zero
            z_rope_upper_on_plot = (options.epsilon_orig - m_fd) / s_fd;
            z_rope_lower_on_plot = (0 - m_fd) / s_fd;
            
            actual_z_rope_lower_on_plot = min(z_rope_lower_on_plot, z_rope_upper_on_plot);
            actual_z_rope_upper_on_plot = max(z_rope_lower_on_plot, z_rope_upper_on_plot);

            % On each subplot (gca):
            y_lim = get(gca, 'YLim');
            fill([actual_z_rope_lower_on_plot, actual_z_rope_upper_on_plot, actual_z_rope_upper_on_plot, actual_z_rope_lower_on_plot],...
                [y_lim(1), y_lim(1), y_lim(2), y_lim(2)], ...
                'm', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'DisplayName', 'ROPE');
%             xline(actual_z_rope_lower_on_plot, 'm--'); xline(actual_z_rope_upper_on_plot, 'm--');
        else
            % Optionally, add text to plot indicating s_fd is ~0 and how to interpret original ROPE
%             text(pos_x, pos_y, sprintf('Original data std is ~0. Values ~%.2f. Original ROPE is +/-%.2f ms', m_fd, options.epsilon_orig));
        end
    elseif strcmp(options.test_type, 'directional_negative')
        if s_fd > sqrt(eps)  && ~isnan(options.epsilon_orig) % Only draw ROPE on z-plot if s_fd is not effectively zero
            z_rope_upper_on_plot = (0 - m_fd) / s_fd;
            z_rope_lower_on_plot = (options.epsilon_orig - m_fd) / s_fd;
            
            actual_z_rope_lower_on_plot = min(z_rope_lower_on_plot, z_rope_upper_on_plot);
            actual_z_rope_upper_on_plot = max(z_rope_lower_on_plot, z_rope_upper_on_plot);

            % On each subplot (gca):
            y_lim = get(gca, 'YLim');
            fill([actual_z_rope_lower_on_plot, actual_z_rope_upper_on_plot, actual_z_rope_upper_on_plot, actual_z_rope_lower_on_plot],...
                [y_lim(1), y_lim(1), y_lim(2), y_lim(2)], ...
                'm', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'DisplayName', 'ROPE');
%             xline(actual_z_rope_lower_on_plot, 'm--'); xline(actual_z_rope_upper_on_plot, 'm--');
        else
            % Optionally, add text to plot indicating s_fd is ~0 and how to interpret original ROPE
%             text(pos_x, pos_y, sprintf('Original data std is ~0. Values ~%.2f. Original ROPE is +/-%.2f ms', m_fd, options.epsilon_orig));
        end
        
    else
        xline(0, 'k--');
    end
    xlim([min(x_hist_range) max(x_hist_range)]);
    xlabel('Normalized Latency Difference (z-score)');
    ylabel('Density');
    legend('show', 'Location', 'best');
    hold off;
    
    sgtitle(['BF (', options.test_type, '): ', num2str(BF)]);
end

% odds_posterior_empirical = posteriorsignedlikelyhood_effect_empirical / (1 - posteriorsignedlikelyhood_effect_empirical);
% % Handle potential division by zero if probability is 1 or 0
% if posteriorsignedlikelyhood_effect_empirical == 1
%     odds_posterior_empirical = Inf;
% elseif posteriorsignedlikelyhood_effect_empirical == 0
%     odds_posterior_empirical = 0; % or handle as very small number / error
% end
% 
% 
% odds_prior_empirical = signedlikelyhood_null_empirical / (1 - signedlikelyhood_null_empirical);
% if signedlikelyhood_null_empirical == 1
%     odds_prior_empirical = Inf;
% elseif signedlikelyhood_null_empirical == 0
%     odds_prior_empirical = 0; % or handle
% end
% 
% % Calculate BF
% if odds_prior_empirical == 0 && odds_posterior_empirical == 0
%     BF = 1; % Or undefined/NaN, as both models predict zero probability for positive values
% elseif odds_prior_empirical == 0 && odds_posterior_empirical > 0
%     BF = Inf; % Strong evidence for effect model
% else
%     BF = odds_posterior_empirical / odds_prior_empirical;
% end
