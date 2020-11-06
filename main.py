from view import View
from core import Core


if __name__ == "__main__":


    #core =  model.PropSpecReaderModel()
    #ctrl = controller.controller(core)
    core = Core()
    view = View(core)
    view.run()