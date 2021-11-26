function iMfield = inv_field(Mfield)
iMfield = -Mfield;
for i = size(Mfield,4)
    iMfield(:,:,:,i,:) = -imwarp4(squeeze(Mfield(:,:,:,i,:)),-Mfield);
end
end