function ds = calculate_ds(locMap,triMap)

ds = Inf;

for i = 1:1:length(triMap)
    v1 = triMap(i,1);
    v2 = triMap(i,2);
    v3 = triMap(i,3);
    
    d1 = (locMap(v1,:)-locMap(v2,:))*(locMap(v1,:)-locMap(v2,:))';
    d2 = (locMap(v1,:)-locMap(v3,:))*(locMap(v1,:)-locMap(v3,:))';
    d3 = (locMap(v3,:)-locMap(v2,:))*(locMap(v3,:)-locMap(v2,:))';
    
    ds = min([d1, d2, d3, ds]);
end;