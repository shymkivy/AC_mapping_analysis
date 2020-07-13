function f_play_stack(stack, fr, loop_num, title_tag)
if ~exist('loop_num', 'var') || isempty(loop_num)
    loop_num = 1;
end
if ~exist('title_tag', 'var') || isempty(title_tag)
    title_tag = '';
end


num_frames = size(stack,3);

figure;
for n_fr = 1:num_frames*loop_num
    n_im_fr = rem(n_fr-1, num_frames)+1;
    imagesc(stack(:,:,n_im_fr)); axis equal tight;
    title(sprintf('%s Frame %d/%d', title_tag, n_im_fr, num_frames))
    drawnow;
    pause(1/fr);
end

end