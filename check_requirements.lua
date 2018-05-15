--- Workaround to emulate try-catch blocks.
-- @see https://www.lua.org/wshop06/Belmonte.pdf
-- @param f function block to run
-- @param catch_f exception handler
function try(f, catch_f)
    local status, exception = pcall(f)
    if not status then
        catch_f(exception)
    end
end

try(function()
    require('torch')
    require('cutorch')
    require('nn')
    require('cunn')
    require('nnx')
    require('cunnx')

    require('optim')

    require('json')
    require('hdf5')
    require('lfs')
    require('os')
    require('io')

    vol_nnup = nn.VolumetricUpSamplingNearest(2, 3, 4)

    print('[Requirements] Seems that all requirements are met. However, this does not guarantee that everything will work.')
    print('[Requirements] Please consult the documentation for details on the requirements and further steps.')
end, function(e)
    print('[Requirements] Seems like one or more requirements are not met, see error mesage below:')
    print(e)
end)